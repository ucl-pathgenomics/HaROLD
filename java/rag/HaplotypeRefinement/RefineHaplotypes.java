import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.util.CombinatoricsUtils;
import pal.alignment.Alignment;
import pal.alignment.AlignmentReaders;
import pal.datatype.Nucleotides;
import picocli.CommandLine;

/**
 * <P>New version of Harold <br>
 * Modified 2 May 2019 
 * 
 * <P>The algorithm takes a set of reads for different samples and calculates the
 * optimal probability of the haplotypes for each of the samples as well as the
 * optimal sequences of the haplotypes
 * <P>We are going to adjust the haplotype frequences and sequences, use these values to assign reads,
 * and then calculate log likelihood based on these assigned reads
 * <P>The way the algorithm works: <br>
 * <blockquote>
 * Initialisation
 * <ul>
 * <li>Reads in parameter values, reference sequence in Fasta format, and BAM files</li>
 * <li>Calculates which sites are variable, which reads contain variable sites, and stores these in different lists</li>
 * <li>Assign reads probabilistically and randomly to haplotypes in randomInitialiseModel by calling Read.initialiseProbHaplo</li>
 * <li>Based on these assignments generate initial values of other parameters</li>
 * </ul>
 * Iteratively
 * <ul>
 * <li>We have current haplotype probabilities and (probabilistic) haplotype sequences</li>
 * <li>With these values we also have the probability that each read is assigned to one of the haplotypes</li>
 * <li>For each sample, we adjust the haplotype probabilities to maximise the log likelihood, where we fix the 
 * haplotype probabilities at the other samples at their current values. This means that when we optimise one sample's
 * haplotype frequencies we 'undo' this optimisation when optimising the haplotype frequencies for the next sample</li>
 * <li>We then, for each variable position, we optimise the haplotype sequences, given the previously optimised haplotype
 * frequencies. Again, we use the unoptimised haplotype sequences at the other positions</li>
 * <li>These optimisations are performed with a term that initially favours non-specific haplotype sequences (e.g. probability of
 * each observed base approximately equal) but then favours specific haplotype sequences (e.g. probability of bases = {0,1})</li>
 * </ul>
 * </blockquote>
 * 

 * @author rgoldst
 */
public class RefineHaplotypes {
    
    int nSites;                                                         // Number of sites
    int nSamples;                                                       // Number of samples
    
    int nRegions;
    
    private boolean[][] basePresent;                                    // Which bases present [iSite][4]
    private int[] nBasePresent;                                         // Number of bases present at site
    private int[][][] siteCount;                                        // [iSample][iSite][strand][5]
    private int[] consensusSeq;
    private ArrayList<Integer> activeSites = new ArrayList<>();         // List of all active (variable) sites
    
    Alignment alignment;
    Alignment refAlignment;
    
    ArrayList<Read> readList = new ArrayList<>();                       // List of reads for each sample
        
    static double errorRate = 0.002;                                    // Expected error rate
    static byte minBaseQual = (byte) 30; // 30;                         // Minimum BAM quality score
    static int minMappingQual = 10;  // 10;                             // Minimum BAM mapping quality
    double minHapFreq = 0.00;
    double minEValue = 1.0E-3;
    double minReadDepth = -1.;
    double minReads = 20.;
    int nOptBaseMaxIter = 1000;
    
    double nParamsPerHaplo = 1.;
    
    boolean iterate = false;
    boolean printRef = false;
    boolean printIntermediate = false;
    int maxIter = -1;
    int maxHaplo = 10;
    int maxRecombine = 1000;
    boolean expand = false;
    
    String outputFileNameTag;
    String sampleTag;
    static PrintWriter alignmentWriter = null;                             // PrintWriter for haplotypes
    static PrintWriter logWriter = null;                                // PrintWriter for log files
    
    double smallValue = 0.;
    double bigValue = 1.;
    char[] dna = {'A', 'C', 'G', 'T', '-', 'X', 'a', 'c', 'g', 't', 'x'}; 

    static Random random;
    MersenneTwister twister;
    static NormalDistribution nd = new NormalDistribution(0., 10.0);
    Options options = new Options();

    
    /**
     * 
     * @param args Parameters and filenames
     */
    public static void main(String[] args) {
        RefineHaplotypes m = new RefineHaplotypes();
        m.run(args);
    }
    
    private RefineHaplotypes() {
        fillHash();
    }

    private void run(String[] args) {
        initialiseOptions(args);
        Parameters startParameters = readData(args);                    // Read in data
        process(startParameters);                                       // Perform refinement
        finish();                                                       // Finish up, close files
    }

    
    void initialiseOptions(String[] args) {
        CommandLine cmd = new CommandLine(options);
        try {
            cmd.parse(args);
            if (cmd.isUsageHelpRequested()) {
                cmd.usage(System.err);
            } else if (cmd.isVersionHelpRequested()) {
                cmd.printVersionHelp(System.err);
            }
        } catch (CommandLine.ParameterException ex) {
            System.err.println(ex.getMessage());
            System.exit(1);
            logWriter.println(ex.getMessage());
            ex.getCommandLine().usage(System.err);
        } catch (Exception ex) {
            throw new CommandLine.ExecutionException(cmd, "Error", ex);
        }
    }
             
    /**
     * Refine the haplotypes in an interative manner
     * 
     * @param startParameters Initial starting configuration based on input data
     */
    public void process(Parameters startParameters) {
        // Optimise parameters for current set of haplotypes
        System.out.println("\nINITIAL OPTIMISATION");
        Parameters bestParameters = optimiseBases(startParameters, 0, 0.0);
        printStuff(bestParameters,0, printIntermediate); 
                
        // Attempt various refinement operations
        // Start with current set of parameters
        Parameters nextBest = bestParameters;  
        // If we are considering recombination events, try random recombinations and see if there is any improvement
        if (maxRecombine > 0) {
            for (int iTry = 0; iTry < maxRecombine; iTry++) {
                // Require at least two haplotypes with non-negligible frequencies
                if (nextBest.nActive > 1) {
                    // Generate state after recombination
                    Parameters recombineParameters = nextBest.recombineHaplo(nSites, activeSites);
                    // Optimise the parameters for this new set of recombined haplotypes
                    recombineParameters = optimiseBases(recombineParameters, 0,
                                    nextBest.adjLogLikelihood);
                    // If penalised log likelihood improves, accept recombination
                    if (recombineParameters.adjLogLikelihood > nextBest.adjLogLikelihood) {
                        nextBest = recombineParameters;
                    }
                    System.out.println("RECOMBINE\t" + iTry + "\t" + recombineParameters.adjLogLikelihood);
                    logWriter.println("RECOMBINE\t" +  iTry + "\t" + recombineParameters.adjLogLikelihood);
                }
            }
            // If improvement, make recombined set into current set
            if (nextBest.adjLogLikelihood > bestParameters.adjLogLikelihood) {
                bestParameters = nextBest;
            }
            System.out.println("INITIAL RECOMBINE\t" + bestParameters.totalLogLikelihood + "\t" 
                    + bestParameters.totPenalty + "\t" + bestParameters.adjLogLikelihood);
            logWriter.println("INITIAL RECOMBINE\t" + bestParameters.totalLogLikelihood + "\t" 
                    + bestParameters.totPenalty + "\t" + bestParameters.adjLogLikelihood);
        }
              
        // If iterate == TRUE, loop around to refine in iterative manner
        if (iterate) {    
            nextBest = bestParameters;
            // Loop over iterations
            for (int bigloop = 1; bigloop < maxIter; bigloop++) {
                
                // See if we can improve model by compressing two haplotypes into one
                // Reduces number of haplotypes by one
                if (bestParameters.nActive > 1) {
                    // Loop over all pairs of haplotypes
                    for (int iHaplo = 0; iHaplo < bestParameters.nHaplo-1; iHaplo++) {
                        for (int jHaplo = iHaplo+1; jHaplo < bestParameters.nHaplo; jHaplo++) {
                            // Make sure that both have non-neglible probabilities
                            if(bestParameters.hapFreq[iHaplo] > 1.0E-6 && 
                                    bestParameters.hapFreq[jHaplo] > 1.0E-6) {
                                // Generate new set of haplotypes; compressing two haplotypes reduces number
                                // of haplotypes by 1
                                Parameters condenseParameters 
                                        = bestParameters.compressHaplo(iHaplo, jHaplo, nSites);
                                // Optimise bases for new set of compressed haplotypes
                                condenseParameters = optimiseBases(condenseParameters, bigloop,
                                    nextBest.adjLogLikelihood);
                                // If penalised log likeklihood improvies, update nextBest to compressed set
                                if (condenseParameters.adjLogLikelihood > nextBest.adjLogLikelihood) {
                                    nextBest = condenseParameters;
                                }
                                System.out.println("CONTRACT\t" + bigloop + "\t" + condenseParameters.adjLogLikelihood);
                                logWriter.println("CONTRACT\t" + bigloop + "\t" + condenseParameters.adjLogLikelihood);
                            }
                        }
                    }
                }
                
                // See if we can improve log likelihood by dividing haplotype into two haplotypes
                // Increases number of haplotypes by one
                if (expand && bestParameters.nHaplo < maxHaplo) {              
                    for (int iHaplo = 0; iHaplo < bestParameters.nHaplo; iHaplo++) {
                        if (bestParameters.hapFreq[iHaplo] > 1.0E-6) {
                            Parameters expandParameters = bestParameters.expandHaplo(iHaplo, nSites);
                            expandParameters = optimiseBases(expandParameters, bigloop,
                                    nextBest.adjLogLikelihood);
                            System.out.println("EXPAND\t" + iHaplo + "\t" + expandParameters.adjLogLikelihood);
                            if (expandParameters.adjLogLikelihood > nextBest.adjLogLikelihood) {
                                nextBest = expandParameters;
                            }
                        }
                    }  
                    
                }

                // If on first loop or there has been an improvement in penalised log likelihood
                // again try random recombinations
                if (bigloop == 1 || nextBest.adjLogLikelihood > bestParameters.adjLogLikelihood) {
                    for (int iTry = 0; iTry < maxRecombine; iTry++) {
                        if (nextBest.nActive > 1) {
                            Parameters recombineParameters = nextBest.recombineHaplo(nSites, activeSites);
                            recombineParameters = optimiseBases(recombineParameters, bigloop,
                                    nextBest.adjLogLikelihood);
                            if (recombineParameters.adjLogLikelihood > nextBest.adjLogLikelihood) {
                                nextBest = recombineParameters;
                            }
                            System.out.println("RECOMBINE\t" + bigloop + "\t" + iTry + "\t" + recombineParameters.adjLogLikelihood);
                            logWriter.println("RECOMBINE\t" + bigloop + "\t" + iTry + "\t" + recombineParameters.adjLogLikelihood);
                        }
                    }
                }
                                    
                System.out.println("\n********** BEST " + bigloop + " " + nextBest.adjLogLikelihood + " *************");
                printStuff(nextBest, bigloop, printIntermediate);
                
                // If there is an improvement, replace bestParameters with best new set of parameters
                if (nextBest.adjLogLikelihood - bestParameters.adjLogLikelihood < 1.0 || nextBest.nActive == 10) {
                    if (nextBest.adjLogLikelihood > bestParameters.adjLogLikelihood) {
                        bestParameters = nextBest;
                    }
                    break;
                } else {
                    bestParameters = nextBest;
                }
                
            }           
        }

        // Output hamming distance between pairs of new haplotypes
        int[][] iMaxBase = bestParameters.iMaxBase;
        for (int iHaplo = 0; iHaplo < bestParameters.nHaplo-1; iHaplo++) {
            for (int jHaplo = iHaplo+1; jHaplo < bestParameters.nHaplo; jHaplo++) {
                if (bestParameters.hapFreq[iHaplo] > 1.0E-6 && bestParameters.hapFreq[jHaplo] > 1.0E-6) {
                    double[] hamming = new double[2];
                    for (int iSite : activeSites) {
                        if(iMaxBase[iHaplo][iSite] < 4 && iMaxBase[jHaplo][iSite]< 4) {
                            hamming[1]++;
                            if (iMaxBase[iHaplo][iSite] != iMaxBase[jHaplo][iSite]) {
                                hamming[0]++;
                            }
                        }
                    }  
                    System.out.println(iHaplo + "\t" + jHaplo + "\t" + (hamming[0]/(1.0E-10+hamming[1])));
                }
            }
        }
        
        logWriter.println("FINISHED");
        if (!printIntermediate) {
            printStuff(nextBest, -1, true);
        }
        
        
        
        logWriter.println("\nHaplotype frequencies");
        for (int iHaplo = 0; iHaplo < bestParameters.nHaplo; iHaplo++) {
            if (bestParameters.hapFreq[iHaplo] > 1.0E-10) {
                logWriter.format("%d\t%.6f\n", iHaplo, bestParameters.hapFreq[iHaplo]);
            }
        }
        logWriter.println("Active sites");
        for (int iSite : activeSites) {
            logWriter.print(iSite);
            for (int iHaplo = 0; iHaplo < bestParameters.nHaplo; iHaplo++) {
                if (bestParameters.hapFreq[iHaplo] > 1.0E-10) {
                    logWriter.format("\t%.6f,%.6f,%.6f,%.6f", bestParameters.probBase[iHaplo][iSite][0],
                            bestParameters.probBase[iHaplo][iSite][1],
                            bestParameters.probBase[iHaplo][iSite][2],
                            bestParameters.probBase[iHaplo][iSite][3]);
                }
            }
            logWriter.println();
        }
        if (printRef) {
            alignmentWriter.println(">" + refAlignment.getIdentifier(0).toString());
            for (int iSite = 1; iSite < nSites; iSite++) {
                alignmentWriter.print(refAlignment.getData(0, iSite-1));
            }
            alignmentWriter.println(); 
        }
    }

    
    
    Parameters optimiseBases(Parameters parameters, int bigLoop, double valToBeat) {
        Parameters nextParameters = assignToHaplotypes(parameters, 0, bigLoop);
        Parameters bestParameters = nextParameters;
        double[] values = new double[nOptBaseMaxIter];
        double bestValue = nextParameters.adjLogLikelihood;
        values[0] = bestValue;
        for (int iIter = 1; iIter < nOptBaseMaxIter; iIter++) {
            nextParameters = assignToHaplotypes(nextParameters, iIter, bigLoop);
            values[iIter] = Math.max(nextParameters.adjLogLikelihood, bestValue);

            if (iIter > 5 && Math.abs(nextParameters.adjLogLikelihood - bestValue) < 1.0) {
                bestParameters = nextParameters;
                break;
            } else if (nextParameters.adjLogLikelihood > bestValue) {
                bestParameters = nextParameters;
                bestValue = bestParameters.adjLogLikelihood;
            }
            if (Math.abs(valToBeat) > 1.0E-9 && iIter > 5) {
                double timeTillBest = (valToBeat-bestValue)/(0.2*(values[iIter]-values[iIter-5]));
                if (timeTillBest > 1000) {
                    break;
                }
            }
        }
        return bestParameters;
    }
    
   
    
    /**
     * 
     * @param parameters
     * @param iTry
     * @param bigLoop
     * @return 
     */
    public Parameters assignToHaplotypes(Parameters parameters, int iTry, int bigLoop) {
        double[][][] probBase = parameters.probBase;                        // Current probability of bases
        double[] hapFreq = parameters.hapFreq;                              // Current haplotype frequencies
        int nHaplo = parameters.nHaplo;                                     // Current number of haplotypes

        //Go over reads assigning to haplotype
        // Calculate log likelihood of each read assigned to each haplotype
        HashMap<Read, double[]> logLikelihoodHash = new HashMap<>();
        for (Read read : readList) {
            double[] logLikelihood = new double[nHaplo];
            HashMap<Integer, Integer> sigSiteHash = read.getSigSiteHash();  // Get list of bases for significant sites
            for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                if (hapFreq[iHaplo] > 1.0E-10) {
                    for (int iSite : sigSiteHash.keySet()) {
                        logLikelihood[iHaplo] 
                                += Math.log(probBase[iHaplo][iSite][sigSiteHash.get(iSite)]);
                    }
                }
            }            
            logLikelihoodHash.put(read, logLikelihood);                     // Add to hash for log likelihood of read given each haplotype
        }
        // Iterate optimising haplotype frequencies
        double[] totalReads = new double[nHaplo];                           // Reads assigned to each haplotype    
        double previousLogLikelihood = -1.0E20;
        boolean finished = false;
        int iFreqIter = 0;
        double[] newHapFreq = Arrays.copyOf(hapFreq, nHaplo);               // Start out with current haplotype frequencies
        
        while (!finished) {  
            double newTotalLogLikelihood = 0.0;
            totalReads = new double[nHaplo];
            for (Read read : readList) {
                double[] logLikelihood = logLikelihoodHash.get(read);         
                // Calculate log probability read is assigned to each haplotype
                double[] logProbHaplo = new double[nHaplo];
                double maxTerm = -1.0E100;
                double iMax = -1;
                // Match to various haplotypes
                for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                    logProbHaplo[iHaplo] = -1.0E10;
                    if (newHapFreq[iHaplo] > 1.0E-10) {
                        logProbHaplo[iHaplo] = Math.log(newHapFreq[iHaplo])
                                + logLikelihood[iHaplo];
                        if (logProbHaplo[iHaplo] > maxTerm) {
                            maxTerm = logProbHaplo[iHaplo];
                            iMax = iHaplo;
                        }
                    }
                }
                double summ = 0.0;
                double[] probHaplo = new double[nHaplo];
                double[] deltaLog = new double[nHaplo];
                double sumProb = 0.0;
                
                // Account for sub-optimal assignments
                for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                    deltaLog[iHaplo] = logProbHaplo[iHaplo] - maxTerm;
                    if (deltaLog[iHaplo] > -10.0) {
                        probHaplo[iHaplo] 
                                = Math.exp(deltaLog[iHaplo]);
                    }
                    summ += probHaplo[iHaplo];
                }
                
                // Calculate log likelihood of assignment to some haplotype
                newTotalLogLikelihood += read.getNCopies() * (maxTerm + Math.log(summ));
                for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                    probHaplo[iHaplo] /= summ;
                    totalReads[iHaplo] += read.getNCopies() * probHaplo[iHaplo];
                }      
            }   // End loop over reads
 
            // Calculate total number of reads assigned to some haplotype
            double tot = 0.0;
            for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                if ( (minReadDepth < 0.0 || totalReads[iHaplo] > (minReadDepth * nSites/250.0)) 
                        && (minReads < 0.0 || totalReads[iHaplo] > minReads) ) {
                    tot += totalReads[iHaplo];
                }
            }     
            // Reconstruct new haplotype frequencies
            for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                newHapFreq[iHaplo] = 0.0;
                // See if site is imformative
                if ( (minReadDepth < 0.0 || totalReads[iHaplo] > (minReadDepth * nSites/250.0)) 
                        && (minReads < 0.0 || totalReads[iHaplo] > minReads) ) {
                    newHapFreq[iHaplo] = totalReads[iHaplo]/tot;
                }                                           // loop over check for suff reads
            }                                               // loop over haplotypes                                            // loop over print haplotypes
            // Check if iteration to adjust frequencies has converged or sufficent attempts tried
            iFreqIter++;
            if (iFreqIter > 100 || Math.abs(newTotalLogLikelihood - previousLogLikelihood) < 1.0) {
                finished = true;
            } else {
                previousLogLikelihood = newTotalLogLikelihood;
            }
        } // End loop over optimising haplotype frequencies
        hapFreq = newHapFreq;
        
        double totalLogLikelihood = 0.0;
        double[][][] reconstructSequence = new double[nHaplo][nSites][4];   // Reconstruction of haplotype sequence    
        totalReads = new double[nHaplo];
        for (Read read : readList) {
            double[] logLikelihood = logLikelihoodHash.get(read);          
            double[] logProbHaplo = new double[nHaplo];
            double maxTerm = -1.0E100;
            double iMax = -1;

            // Match to various haplotypes
            for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                logProbHaplo[iHaplo] = -1.0E10;
                if (hapFreq[iHaplo] > 1.0E-10) {
                    logProbHaplo[iHaplo] = Math.log(hapFreq[iHaplo])
                            + logLikelihood[iHaplo];
                    if (logProbHaplo[iHaplo] > maxTerm) {
                        maxTerm = logProbHaplo[iHaplo];
                        iMax = iHaplo;
                    }
                }
            }
            double summ = 0.0;
            double[] probHaplo = new double[nHaplo];

            double[] deltaLog = new double[nHaplo];
            double sumProb = 0.0;
            for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                deltaLog[iHaplo] = logProbHaplo[iHaplo] - maxTerm;
                if (deltaLog[iHaplo] > -10.0) {
                    probHaplo[iHaplo] 
                            = Math.exp(deltaLog[iHaplo]);
                }
                summ += probHaplo[iHaplo];
            }

            // 
            totalLogLikelihood += read.getNCopies() * (maxTerm + Math.log(summ));
            for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                probHaplo[iHaplo] /= summ;
                totalReads[iHaplo] += read.getNCopies() * probHaplo[iHaplo];
            }

            for (int iSite : read.getSigSiteHash().keySet()) {
                int iSeq = read.getSigSiteHash().get(iSite);
                for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                    reconstructSequence[iHaplo][iSite][iSeq] 
                            += probHaplo[iHaplo] * read.getNCopies();
                }
            }
        }  

        double tot = 0.0;
        for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
            if ( (minReadDepth < 0.0 || totalReads[iHaplo] > (minReadDepth * nSites/250.0)) 
                    && (minReads < 0.0 || totalReads[iHaplo] > minReads) ) {
                tot += totalReads[iHaplo];
            }
        }
        
        // Reconstruct new haplotypes
        int nTotalAssigned = 0;
        int[] nAssigned = new int[nHaplo];
        ArrayList<Integer>[] polySites = new ArrayList[nHaplo];             // Array of polymorphic sites
        int[] nPoly = new int[nHaplo];
        double[][][] newProbBase = new double[nHaplo][nSites][4];     
        int[][] iMaxBase = new int[nHaplo][nSites];                         // Sequence of haplotype
        
        for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
            polySites[iHaplo] = new ArrayList<>();
            newHapFreq[iHaplo] = 0.0;
            Arrays.fill(iMaxBase[iHaplo], 4);
            for (int iSite : activeSites) {
                Arrays.fill(newProbBase[iHaplo][iSite], 0.25);
            }
            // See if site is imformative
            if ( (minReadDepth < 0.0 || totalReads[iHaplo] > (minReadDepth * nSites/250.0)) 
                    && (minReads < 0.0 || totalReads[iHaplo] > minReads) ) {
                newHapFreq[iHaplo] = totalReads[iHaplo]/tot;
//                    alignmentWriter.println(">" + sampleTag + "_H" + (iHaplo+1) + "_" + iTry);
                for (int iSite : activeSites) {
                    double summ = reconstructSequence[iHaplo][iSite][0]
                            + reconstructSequence[iHaplo][iSite][1]
                            + reconstructSequence[iHaplo][iSite][2]
                            + reconstructSequence[iHaplo][iSite][3];
                    if (summ > 0.99999) {
                        double[] fracBase = new double[5];
                        for (int iBase = 0; iBase < 4; iBase++) {
                            double fracObsBase = reconstructSequence[iHaplo][iSite][iBase]/summ;
                            if (fracObsBase > errorRate) {
                                fracBase[iBase] = (fracObsBase-errorRate)/(1.0 - 4.0 * errorRate);
                                fracBase[4] += fracBase[iBase];
                            }
                        }
                        for (int iBase = 0; iBase < 4; iBase++) {
                            newProbBase[iHaplo][iSite][iBase] = 
                                    ((1.0 - 3.0 * errorRate) * (fracBase[iBase]/fracBase[4]) 
                                         + errorRate * (1.0 - (fracBase[iBase]/fracBase[4])));
                            if ((reconstructSequence[iHaplo][iSite][iBase]/summ) > 0.6) {
                                nTotalAssigned++;
                                nAssigned[iHaplo]++;
                                iMaxBase[iHaplo][iSite] = iBase;
                            }
                        } 
                        if (iMaxBase[iHaplo][iSite] == 4) {
                            nPoly[iHaplo]++;
                            polySites[iHaplo].add(iSite);
                        }
                    }
                }
            }                                           // loop over check for suff reads
        }                                               // loop over haplotypes                                            // loop over print haplotypes

        double totPenalty = nHaplo * nParamsPerHaplo;
          
        System.out.println("\n" + iTry + "\t" +  nTotalAssigned + "\t" + 
                totalLogLikelihood + "\t" + totPenalty + "\t" + (totalLogLikelihood - totPenalty));
        for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
            if (hapFreq[iHaplo] > 0.001) {
            System.out.println("\t" + iHaplo + "\t" + hapFreq[iHaplo] + "\t"
                    + nAssigned[iHaplo] + "\t" + nPoly[iHaplo] + "\t"
                    + (nAssigned[iHaplo] + nPoly[iHaplo]) + "\t" + totalReads[iHaplo]);
            }
        }
        Parameters newParameters = new Parameters(totalLogLikelihood, totPenalty, newHapFreq, newProbBase, 
                iMaxBase, polySites);
        return newParameters;
    }

    void printStuff(Parameters inputParameters, int iTry, boolean printAlignment) {
        int[][] iMaxBase = inputParameters.iMaxBase;
        if (iTry == 0) {
            System.out.print("\nINITIAL");
            logWriter.print("\nINITIAL");
        } else {
            System.out.print("\nBig Iteration " + iTry);
            logWriter.print("\nBig Iteration " + iTry);
        }
        System.out.println("\tLL: " + inputParameters.totalLogLikelihood + "\tPENALTY: " 
                + inputParameters.totPenalty + "\tPENALISED LL: " + inputParameters.adjLogLikelihood);
        logWriter.println("\tLL: " + inputParameters.totalLogLikelihood + "\tPENALTY: " 
                + inputParameters.totPenalty + "\tPENALISED LL: " + inputParameters.adjLogLikelihood);

        double[] hapFreq = inputParameters.hapFreq;
        for (int iHaplo = 0; iHaplo < inputParameters.nHaplo; iHaplo++) {
            if (hapFreq[iHaplo] > 1.0E-6) {
                System.out.println("\t" + iHaplo + "\t" + hapFreq[iHaplo]);
                logWriter.println("\t" + iHaplo + "\t" + hapFreq[iHaplo]);
                if (printAlignment) {
                    if (iTry >= 0) {
                        alignmentWriter.println(">" + sampleTag + "_H" + iHaplo + "_" + iTry);
                    } else {
                        alignmentWriter.println(">" + sampleTag + "_H" + iHaplo);
                    }
                    for (int iSite = 1; iSite < nSites; iSite++) {
                        if (activeSites.contains(iSite)) {
                            alignmentWriter.print(dna[iMaxBase[iHaplo][iSite]]);
                        } else {
                            alignmentWriter.print(dna[consensusSeq[iSite]]);
                        }
                    }
                    alignmentWriter.println();
                }
            }
        }
        alignmentWriter.flush();
        System.out.println();
        logWriter.println();
        
    }
    

    /**
     * Read in parameter values, reference sequence, and BAM files
     * 
     * @param args Arguments in String format, containing initial number of
     * haplotypes, number of sites, and name of input file
     */
    public Parameters readData(String[] args) {

        if (options.randomSeed < 0.0) {
            random = new Random();
            nd = new NormalDistribution(0., 10.0);
        } else {
            random = new Random(options.randomSeed);
            twister = new MersenneTwister(options.randomSeed);
            nd = new NormalDistribution(twister, 0., 10.0);
        }
        int nHaplo = 0;
        double[] hapFreq = null;
        double[][][] probBase = null;
        this.outputFileNameTag = options.tag;
        this.minReadDepth = options.minReadDepth;
        this.minReads = options.minReads;
        this.errorRate = options.errorRate; 
        this.iterate = options.iterate;
        this.maxIter = options.maxIter;
        this.maxHaplo = options.maxHaplo;
        this.maxRecombine = options.maxRecombine;
        this.expand = options.expand;
        this.printRef = options.printRef;
        this.printIntermediate = options.printIntermediate;
        openFiles();
        logWriter.print("Command arguments: ");
        System.out.print("Command arguments: ");
        for (int iArg = 0; iArg < args.length; iArg++) {
            logWriter.print(" " + args[iArg]);
            System.out.print(" " + args[iArg]);  
        }


        try {
            // Read in alignment of haplotypes
            Nucleotides dataType = new Nucleotides();
            System.out.println(options.hapAlignmentFile.getAbsolutePath());
            FileReader in = new FileReader(options.hapAlignmentFile);
            alignment = AlignmentReaders.readFastaSequences(in, dataType);
            in.close();
            nHaplo = alignment.getSequenceCount();
            
            // Read in reference alignment
            in = new FileReader(options.refSeqFile);
            refAlignment = AlignmentReaders.readFastaSequences(in, dataType);
            in.close();            
            nSites = refAlignment.getSiteCount()+1;
            nRegions = nSites / 1000;
            
            // Set parameters to appropriate dimensions
            siteCount = new int[nSites][3][5];
            hapFreq = new double[nHaplo];
            probBase = new double[nHaplo][nSites][4];
            consensusSeq = new int[nSites];
            Arrays.fill(consensusSeq, 4);
            FileReader file = null;
            BufferedReader buff = null;
            String line = "";
            String[] words = null;

            if (options.hapFreqFile.exists()) {
                // Read in haplotype frequency file
                file = new FileReader(options.hapFreqFile);
                buff = new BufferedReader(file);
                line = buff.readLine();
                words = line.split("\\t");
                for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                    hapFreq[iHaplo] = Double.parseDouble(words[iHaplo+1]);
                }
                buff.close();
                file.close();
            } else {
                System.out.println("Initialising hapfreq with even values");
                 for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                    hapFreq[iHaplo] = 1.0/nHaplo;
                }               
            }
            
            // Read in probabilities of bases for haplotypes
            file = new FileReader(options.baseFreqFile);
            buff = new BufferedReader(file);
            boolean eof = false;
            while (!eof) {
                line = buff.readLine();
                if (line == null) {
                    eof = true;
                } else {
                    words = line.split("\\t");
                    int iSite = Integer.parseInt(words[0]);
                    for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                        String[] freqs = words[iHaplo+4].split(",");
                        for (int iBase = 0; iBase < 4; iBase++) {
                            double rawProb = Double.parseDouble(freqs[iBase]);
                            probBase[iHaplo][iSite][iBase] 
                                    = rawProb * (1.0 - 3.0 * errorRate) + errorRate * (1.0 - rawProb);
                        }
                    }
                }
            }            
            
            // Read in BAM files
                File sampleFile = options.readsFile;
                sampleTag = sampleFile.getName().replaceAll(".*\\/", "").replaceAll(".bam", "").replaceAll(".sam", "").replaceAll(".BAM", "").replaceAll(".SAM", "");
                readFromFile(sampleFile);
        } catch (IOException e) {
            System.out.println("Error: File not found (IO error)");
            System.exit(1);
        }
       
        basePresent = new boolean[nSites][4];                       // What bases are present at each site
        nBasePresent = new int[nSites];

        for (int iSite = 0; iSite <  nSites; iSite++) {
            // Check for uneven distribution amongst strands
            double[] pValues = checkStrandCount(siteCount[iSite]);
            if (Double.isNaN(pValues[1]) || pValues[1]*nSites > minEValue) {
                int iMax = -1;                                  // Index of most common base
                int nMax = 0;                                   // Number of most common base
                int nTot = 0;                                   // Total number of bases
                for (int iBase = 0; iBase < 4; iBase++) {
                    nTot += siteCount[iSite][2][iBase];
                    if (siteCount[iSite][2][iBase] > nMax) {
                        iMax = iBase;
                        nMax = siteCount[iSite][2][iBase];
                    }
                }
                // See which bases are present
                if (nMax >= 1.999) {         // Require at least 2 total
                    basePresent[iSite][iMax] = true;       // Most common base is present
                    consensusSeq[iSite] = iMax;
                    double noise = Math.max(1.999, errorRate * nMax);         // Require at least 2 and over noise for other bases
                    for (int iBase = 0; iBase < 4; iBase++) {
                        // For each base, see if there are sufficient examples
                        if (iBase != iMax && siteCount[iSite][2][iBase] > noise) {
                                basePresent[iSite][iBase] = true;
                        }
                    }
                }
            }
        }      

        for (int iSite = 0; iSite < nSites; iSite++){
            for (int iBase = 0; iBase < 4; iBase++) {
                if (basePresent[iSite][iBase]) {
                    nBasePresent[iSite]++;
                }
            }
            if (nBasePresent[iSite] > 1) {                                     // Variable site
                nParamsPerHaplo += nBasePresent[iSite]-1.0;
                activeSites.add(iSite);
            }
        }
        System.out.println(nSites + " sites, " + activeSites.size() + " variable");
        int nReads = readList.size();
        
        HashMap<String, Read> compressReadListHash = new HashMap<>();
        for (Read read : readList) {
            read.setSignificantSites(activeSites);
            String tag = read.getSigSiteTag();
            if (compressReadListHash.containsKey(tag)) {
                compressReadListHash.get(tag).addCopies();
            } else if (tag.length() > 0) {
                compressReadListHash.put(tag, read);
            }
        }
        readList = new ArrayList(compressReadListHash.values());
        System.out.println(nReads + " reads, " + readList.size() + " relevant");
        
        Parameters parameters = new Parameters(0.0, 0.0, 
                hapFreq, probBase, null, null);
        
        return parameters;
        
        
    }
    
    /**
     * Reads from BAM file, saves in readList, and updates siteCount
     * 
     * @param inputSamOrBamFile BAM file
     * @param iSample Sample number
     * @throws IOException 
     */
    private void readFromFile(File inputSamOrBamFile) throws IOException {
        SamReader reader = SamReaderFactory.makeDefault().open(inputSamOrBamFile);
        for (final SAMRecord samRecord : reader) { 
            if (samRecord.getMappingQuality() >= RefineHaplotypes.minMappingQual) {
                Read newRead = new Read(samRecord);
                readList.add(newRead);
                int[] limits = newRead.getLimits();
                int[] sequence = newRead.getSequence();
                boolean[] siteExists = newRead.getSiteExists();
                for (int iSite = limits[0]; iSite < limits[1]; iSite++) {
                    int lSite = iSite - limits[0];
                    if (siteExists[lSite]) {
                        siteCount[iSite][2][sequence[lSite]]++; 
                        if (newRead.getNegativeStrand()) {
                            siteCount[iSite][1][sequence[lSite]]++; 
                        } else {
                            siteCount[iSite][0][sequence[lSite]]++; 
                        }
                    }
                }
            }
        }        
    }
    
    
    double[] checkStrandCount(int[][] sCount) {
        double[] pValues = new double[2];
        int[] iCount = new int[2];
        iCount[0] = sCount[0][0]
                + sCount[0][1]
                + sCount[0][2]
                + sCount[0][3];
        iCount[1] = sCount[1][0]
                + sCount[1][1]
                + sCount[1][2]
                + sCount[1][3];            
        int total = iCount[0] + iCount[1]; 
        int iMin = 0;
        if (iCount[1] < iCount[0]) {
            iMin = 1;
        }
        BinomialDistribution bd = new BinomialDistribution(total, 0.5);
        pValues[0] = bd.cumulativeProbability(iCount[iMin]);
        pValues[1] = 100.0;

        if (iCount[0] > 2 && iCount[1] > 2) {
            for (int iBase = 0; iBase < 4; iBase++) {
                int sumBase = sCount[0][iBase]+sCount[1][iBase];
                if (sumBase > 0 && sumBase < (total+1)/2) {
                    double realTerm = CombinatoricsUtils.binomialCoefficientDouble(iCount[0], sCount[0][iBase])
                            * CombinatoricsUtils.binomialCoefficientDouble(iCount[1], sCount[1][iBase]);
                    double[] p = new double[2];
                    for (int iTerm = 0; iTerm <= Math.min(sumBase, iCount[0]); iTerm++) {
                        if (sumBase-iTerm <= iCount[1]) {
                            double term 
                                    = CombinatoricsUtils.binomialCoefficientDouble(iCount[0], iTerm)
                                    * CombinatoricsUtils.binomialCoefficientDouble(iCount[1], (sumBase-iTerm));
                            if (term <= realTerm+0.00001) {
                                p[0] += term;
                            } else {
                                p[1] += term;
                            }
                        }
                    }
                    pValues[1] = Math.min(p[0]/(p[0]+p[1]), pValues[1]);
                }
            }
        }
        return pValues;
    }
    
    
    
    
    static private final HashMap<Character, Integer> DNA_CHAR_TO_INT_DNA_HASH = new HashMap<>();
    private void fillHash() {
        DNA_CHAR_TO_INT_DNA_HASH.put('A', 0);
        DNA_CHAR_TO_INT_DNA_HASH.put('C', 1);
        DNA_CHAR_TO_INT_DNA_HASH.put('G', 2);
        DNA_CHAR_TO_INT_DNA_HASH.put('T', 3);
        DNA_CHAR_TO_INT_DNA_HASH.put('-', 4);
        DNA_CHAR_TO_INT_DNA_HASH.put('N', 4);
    }    

    /**
     * Converts DNA in character format to DNA in integer format
     * @param DNA DNA in character format
     * @return DNA in integer format
     */
    public static int getDNA(char DNA) {
        return DNA_CHAR_TO_INT_DNA_HASH.get(DNA);
    }
    
    public void openFiles() {
    try{
        File logFile = new File(outputFileNameTag + ".log");
        if (logFile.exists()){
            logFile.delete();
        }
        logWriter = new PrintWriter(logFile);    

        File alignmentOutputFile = new File(outputFileNameTag + ".fasta");
        if (alignmentOutputFile.exists()){
            alignmentOutputFile.delete();
        }
        alignmentWriter = new PrintWriter(alignmentOutputFile);                         
    } catch (IOException ioe) {
        System.out.println("Error -- " + ioe.toString());
        logWriter.println("Error -- " + ioe.toString());
        System.exit(1);
        }
    }
    
    
        
    public void finish() {
        try{
            logWriter.close();
            alignmentWriter.close();        
        } catch(Exception ex) {
            System.exit(1);
        }
    }
    
}