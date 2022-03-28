/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package cluster_rg;

import flanagan.math.MinimisationFunction;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;


/**
 *
 * @author rgoldst
 */
public class DataSet {
    int nTimePoints = 0;   // Number of time points

    private ArrayList<Site> activeSiteVector = new ArrayList<>();  // List of sites that are actively considered
    private ArrayList<Site> variableSiteVector = new ArrayList<>(); // List of all variable sites
    private ArrayList<Site> reducedSiteVector0 = new ArrayList<>();
    private ArrayList<Site> reducedSiteVector1 = new ArrayList<>();
    private int nHaplo = 3; // Number of haplotypes
    private ArrayList<Assignment> assignmentVector = null;   // Vectir if assignments
    private int[] nAssignDiffBases = null;
    private double[] currentAlphaParams = new double[2];   // alpha0 and alphaE
    private double[][] currentPiHap = null;
    private double[] useFrac = Constants.USE_FRAC;
    private int iIter = 0;
    private double[] priors = null;
    private boolean verbose;
    double[] avgReadDepth = null;

    private int[] nActive = null;
    private double[][] fixedValues = null;
    private int[][] mapping = null;   
    private double[][] initialHapParams = null;
    final int siteCount;
    GammaCalc gammaCalc = null;
    Random random = null;
    int nEvaluations = 0;
    
    private String[] baseString = {"A", "C", "G", "T", " ", "-"};
        
    private int optType = 0;  // 0 for optimising alpha0 and alphaE, 1 for optimising haplotype frequencies
    private int optTimePoint = 0;   // if optType = 1, what timePoint is being optimised
    private int iCount = 0;  // How many iterations of optimiser have been finished
    
    private double currentLogLikelihood = 0.0;

    DataSet(File fileNameFile, File initialFreqFile, int nHaplo, ArrayList<Assignment> assignmentVector,
            int[] nAssignDiffBases, GammaCalc gammaCalc, Random random, boolean verbose) {  // Read in data
        this.nHaplo = nHaplo;
        this.assignmentVector = assignmentVector;
        this.nAssignDiffBases = nAssignDiffBases;
        this.verbose = verbose;
        this.gammaCalc = gammaCalc;
        this.random = random;
        priors = new double[5];    
        ArrayList<Integer> allSiteVector = new ArrayList<>();// List of all sites
        HashMap<Integer, Site> siteHash = new HashMap<Integer, Site>();  // Data of sites labeled by site number
        priors[1] = Math.log(0.9 / (nAssignDiffBases[1]+1.0E-20));
        priors[2] = Math.log(0.07 / (nAssignDiffBases[2]+1.0E-20));
        priors[3] = Math.log(0.02 / (nAssignDiffBases[3]+1.0E-20));
        priors[4] = Math.log(0.01 / (nAssignDiffBases[4]+1.0E-20));
        List<String> fileNameVector;
        List<String> initialFreqVector = new ArrayList<>();
        boolean initialFreqSet = false;

        try {
            System.out.println(fileNameFile);
            Cluster_RG.logWriter.println(fileNameFile);
            fileNameVector = Files.readAllLines(fileNameFile.toPath());
            if (!(initialFreqFile == null)) {
                initialFreqVector = Files.readAllLines(initialFreqFile.toPath());
                initialFreqSet = true;
            }
        }
        catch (IOException e) {
            System.out.println("Error: File not found (IO error)");
            System.out.println(fileNameFile);
            throw new RuntimeException(e);
        }

        nTimePoints = fileNameVector.size();  // Number of timepoints = number of files
        initialiseHapParams(fileNameFile, initialFreqSet, initialFreqVector);  
        
        String[] dataFile = new String[nTimePoints];
        String pathPrefix = Paths.get(fileNameFile.getAbsolutePath()).getParent().toString();
        for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {  // read in data files
            dataFile[iTimePoint] = new File(pathPrefix, fileNameVector.get(iTimePoint)).toString();
            try{
                System.out.println(dataFile[iTimePoint]);
                Cluster_RG.logWriter.println(dataFile[iTimePoint]);
                FileReader file = new FileReader(dataFile[iTimePoint]);
                BufferedReader buff = new BufferedReader(file);
                String line = "";
                boolean eof = false;
                while (!eof) {
                    line = buff.readLine();
                    if (line == null) {
                            eof = true;
                    } else {
                        if (line.contains("Position")) {
                            line = buff.readLine();
                        }
                        int iSite = Integer.parseInt(line.split(",")[0]);
                        if (!siteHash.containsKey(iSite)) {   // list of sites that contain data
                            allSiteVector.add(iSite);
                            Site newSite = new Site(iSite, nTimePoints, nActive, nHaplo, assignmentVector, gammaCalc); // create new site if needed
                            siteHash.put(iSite, newSite);
                        }
                        siteHash.get(iSite).addTimePoint(iTimePoint, line);  // add datapoint to site
                        
                    }
                }
            }     
            catch (IOException e) {
                System.out.println("Error: File not found (IO error)");
                System.exit(1);
            }

        }
        this.siteCount = siteHash.size();
        currentPiHap = new double[nTimePoints][nHaplo];
        avgReadDepth = new double[nTimePoints];

        for (int iSite : allSiteVector) {  // Create activeSiteVector
            Site site = siteHash.get(iSite);
            for(int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                int[] totReads = site.getTotReads();
                avgReadDepth[iTimePoint] += 1.0 * totReads[iTimePoint];
            }
            if (site.isActive()) {    // do simple sums
                activeSiteVector.add(site);
                if (random.nextDouble() < useFrac[0]) {
                    reducedSiteVector0.add(site);
                }
                if (random.nextDouble() < useFrac[1]) {
                    reducedSiteVector1.add(site);
                }
                if (!site.siteConserved) {
                    variableSiteVector.add(site);
                }
            } 
        }
        System.out.println("Average read depth");
        Cluster_RG.logWriter.println("Average read depth");
        for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
            System.out.println(iTimePoint + "\t" + avgReadDepth[iTimePoint]/allSiteVector.size());
            Cluster_RG.logWriter.println(iTimePoint + "\t" + avgReadDepth[iTimePoint]/allSiteVector.size());
        }
        System.out.println();
        Cluster_RG.logWriter.println();
    }
    
    double computeTotalLogLikelihood(double[] alphaParams, double[][] piHap) {
        currentAlphaParams = alphaParams;
        currentPiHap = piHap;
        double alpha0 = alphaParams[0]*(1.0-alphaParams[1])/alphaParams[1];
        double alphaE = (1.0-alphaParams[0])*(1.0-alphaParams[1])/alphaParams[1];
        double lGA0 = this.gammaCalc.logGamma(alpha0);
        double lGA03AE = this.gammaCalc.logGamma(alpha0 + 3.0 * alphaE);

        if (optType == 0 || optType == 2) {
            nEvaluations += nTimePoints;
            double totalLogLikelihood = 0.0;
            for (Site site : activeSiteVector) {
                totalLogLikelihood += 
                        site.assignHaplotypes(alpha0, alphaE, lGA0, lGA03AE, priors, piHap);
            }
            return totalLogLikelihood;
        } else if (optType == 1) {
            nEvaluations++;
            double totalLogLikelihood = 0.0;
            for (Site site : activeSiteVector) {
                totalLogLikelihood += 
                        site.computeSiteTimePointLogLikelihood(optTimePoint, 
                        alpha0, alphaE, lGA0, lGA03AE, priors, piHap[optTimePoint]);
            }
            return totalLogLikelihood;     
        }
        return 0.0;
    } 
    
    void setOptType(int optType, int optTimePoint, int iIter) {
        this.optType = optType;
        this.optTimePoint = optTimePoint;
        this.iIter = iIter;
        if (this.verbose) {
            System.out.println("ggg\t" + iIter + "\t" + optType + "\t" + optTimePoint);
        }
        iCount = 0;
    }
    

    void updatePriors() {
        double[] nPresent = new double[Constants.MAX_BASES+1];
        Arrays.fill(nPresent, 5.0);
        for (Site site : activeSiteVector) {
            nPresent[site.getBestGuessNPresent()]++;
        }
        double summ = nPresent[1]+nPresent[2]+nPresent[3]+nPresent[4];
        nPresent[1] /= summ;
        nPresent[2] /= summ;
        nPresent[3] /= summ;
        nPresent[4] /= summ;
        priors[1] = Math.log(nPresent[1] / (nAssignDiffBases[1]+1.0E-20));
        priors[2] = Math.log(nPresent[2] / (nAssignDiffBases[2]+1.0E-20));
        priors[3] = Math.log(nPresent[3] / (nAssignDiffBases[3]+1.0E-20));
        priors[4] = Math.log(nPresent[4] / (nAssignDiffBases[4]+1.0E-20));
    }
    
    public double value(double[] params) {     
        if (optType == 0) {
            double val = computeTotalLogLikelihood(params, currentPiHap);
            if (iCount % 10 == 0 && this.verbose) {
                System.out.print(Arrays.toString(currentAlphaParams));
                for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                        System.out.print("  " + Arrays.toString(currentPiHap[iTimePoint]));
                }
                System.out.println();
                System.out.print("yyy\t" + Arrays.toString(params));
                System.out.println("\t" + val);
            }
            iCount++;
            return -val;
        } else if (optType == 1) {
            currentPiHap[optTimePoint] = computePiHap(optTimePoint,params);
            double val = computeTotalLogLikelihood(currentAlphaParams, currentPiHap);
            if (this.verbose && iCount % 10 == 0) {
                System.out.print(Arrays.toString(currentAlphaParams));
                for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                        System.out.print("  " + Arrays.toString(currentPiHap[iTimePoint]));
                }
                System.out.println();
                System.out.print("xxx\t" + Arrays.toString(params));
                System.out.println("\t" + val);
            }
            iCount++;
            return -val;
        }
        System.out.println("Error in optimisation");
        System.exit(1);
        return 0.0;
    }
    
    public double value(double singleParam) {  
        double[] params = new double[1];
        params[0] = singleParam;
        return value(params);
    }

    public double value() {
        System.out.println(Arrays.toString(currentAlphaParams));
        for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
            System.out.println(Arrays.toString(currentPiHap[iTimePoint]));
        }
        return computeTotalLogLikelihood(currentAlphaParams, currentPiHap);
    }
    
    void setParams(double[][] currentHapParams, double[] currentAlphaParams) {
        for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {;
            currentPiHap[iTimePoint] = computePiHap(iTimePoint, currentHapParams[iTimePoint]);
        }
        this.currentAlphaParams =  currentAlphaParams;
    }
    
    void printResults(boolean printHaplotypes, boolean printLogLikeli) {
        int nParams = 3 + (nHaplo-1)*nTimePoints;
        System.out.println("\nResults for nHaplotypes = " + nHaplo);
        System.out.println("Number of adjustable parameters: " + nParams);
        System.out.println("Final likelihood: " + currentLogLikelihood);
        System.out.println("Dirichlet parameters for errors: " + currentAlphaParams[0] + "\t" + currentAlphaParams[1]
                + "\t" + (nAssignDiffBases[1] * Math.exp(priors[1])) 
                + "\t" + (nAssignDiffBases[2] * Math.exp(priors[2]))
                + "\t" + (nAssignDiffBases[3] * Math.exp(priors[3]))
                + "\t" + (nAssignDiffBases[4] * Math.exp(priors[4])) );
        System.out.println("Error rate: " + (1.0 - currentAlphaParams[0]));
        System.out.println("StdDev: " + Math.sqrt(currentAlphaParams[1]*currentAlphaParams[0]*(1.0 - currentAlphaParams[0])));
        System.out.println("Haplotype frequencies");
        Cluster_RG.logWriter.println("\nResults for nHaplotypes = " + nHaplo);
        Cluster_RG.logWriter.println("Number of adjustable parameters: " + nParams);
        Cluster_RG.logWriter.println("Final likelihood: " + currentLogLikelihood);
        Cluster_RG.logWriter.println("Dirichlet parameters for errors: " + currentAlphaParams[0] + "\t" + currentAlphaParams[1]
                + "\t" + (nAssignDiffBases[1] * Math.exp(priors[1])) 
                + "\t" + (nAssignDiffBases[2] * Math.exp(priors[2]))
                + "\t" + (nAssignDiffBases[3] * Math.exp(priors[3]))
                + "\t" + (nAssignDiffBases[4] * Math.exp(priors[4])) );
        Cluster_RG.logWriter.println("Error rate: " + (1.0 - currentAlphaParams[0]));
        Cluster_RG.logWriter.println("StdDev: " + Math.sqrt(currentAlphaParams[1]*currentAlphaParams[0]*(1.0 - currentAlphaParams[0])));
        Cluster_RG.logWriter.println("Haplotype frequencies");        
        for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
            System.out.print(iTimePoint);
            Cluster_RG.logWriter.print(iTimePoint);
            for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                Cluster_RG.logWriter.print("\t" + currentPiHap[iTimePoint][iHaplo]);          
            }
            System.out.println();
            Cluster_RG.logWriter.println();
        }
        
        if (printHaplotypes || printLogLikeli) {
            int nSites = 0;
            for (Site site : activeSiteVector) {
                nSites = Math.max(nSites, site.iSite+1);
            }
            int[][] bestBase = new int[nHaplo][nSites];
            double[][] probBestBase = new double[nHaplo][nSites];
            for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                for (int iSite = 0; iSite < nSites; iSite++) {
                    bestBase[iHaplo][iSite] = 5;
                    probBestBase[iHaplo][iSite] = 0.5;
                }
            }

            for (Site site : activeSiteVector) {
                int iSite = site.iSite;
                if (site.siteConserved) {
                    for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                        bestBase[iHaplo][iSite] = site.conservedBase;
                        probBestBase[iHaplo][iSite] = 1.0;
                    }
                } else {
                    double[][] probBase = site.getProbBase();
                    for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                        for (int iBase = 0; iBase < 4; iBase++) {
                            if (probBase[iHaplo][iBase] > probBestBase[iHaplo][iSite]) {
                                probBestBase[iHaplo][iSite] = probBase[iHaplo][iBase];
                                bestBase[iHaplo][iSite] = iBase;
                            }
                        }
                    }
                }
            }
            if (printHaplotypes) {
                for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                    Cluster_RG.hapWriter.println(">Haplo_" + iHaplo);
                    for (int iSite = 1; iSite < nSites; iSite++) {
    //                    for (int iSite = 230000; iSite < 230880; iSite++) {
                        Cluster_RG.hapWriter.print(baseString[bestBase[iHaplo][iSite]]);
                        if (iSite % 80 == 0) {
                            Cluster_RG.hapWriter.println();
                        }
                    }
                    Cluster_RG.hapWriter.println();
                }
            }
            if (printLogLikeli) {
                boolean printStuff = false;
                double alpha0 = currentAlphaParams[0]*(1.0-currentAlphaParams[1])/currentAlphaParams[1];
                double alphaE = (1.0-currentAlphaParams[0])*(1.0-currentAlphaParams[1])/currentAlphaParams[1];
                double lGA0 = this.gammaCalc.logGamma(alpha0);
                double lGA03AE = this.gammaCalc.logGamma(alpha0 + 3.0 * alphaE);
//                for (Site site : activeSiteVector) {
//                    site.computeExpectedLikelihood(alpha0, alphaE, lGA0, lGA03AE, priors);
//                }
                for (Site site : activeSiteVector) {
                    Cluster_RG.logLikeWriter.format("%d\t%.4f\t%.4f", site.iSite, site.totalLogLikelihood, site.expectedLikelihood); 
                    Cluster_RG.logLikeWriter.print("\t" + Arrays.toString(site.strandReads[0][0]) + Arrays.toString(site.strandReads[0][1]));
                    double[][] probBase = site.getProbBase();
                    for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                        Cluster_RG.logLikeWriter.format("\t%.4f,%.4f,%.4f,%.4f", 
                                probBase[iHaplo][0], probBase[iHaplo][1], probBase[iHaplo][2], probBase[iHaplo][3]);
                    }
                    Cluster_RG.logLikeWriter.println();
                }
            }

        }
        
        
        
    }
    

        /**
     *
     * Compute new values of piHap for all time points
     */    
    double[][] computePiHap(double[][] hapParams) {
        double[][] piHap = new double[nTimePoints][nHaplo];
        for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
            piHap[iTimePoint] = computePiHap(iTimePoint, hapParams[iTimePoint]);
        }
        return piHap;
    }
    
    /**
     *
     * compute new values of piHap for single time point
     */    
    double[] computePiHap(int iTimePoint, double[] hapParams) {
        double[] expandedHapParams = Arrays.copyOf(fixedValues[iTimePoint], nHaplo-1);
        for (int iParam = 0; iParam < nActive[iTimePoint]; iParam++) {
            expandedHapParams[mapping[iTimePoint][iParam]] = hapParams[iParam];
        }
        double[] piHap = new double[nHaplo];
        double remaining = 1.0;
        for (int iHaplo = 0; iHaplo < nHaplo-1; iHaplo++) {
            piHap[iHaplo] = remaining * expandedHapParams[iHaplo];
            remaining -= piHap[iHaplo];
        }
        piHap[nHaplo-1] = remaining;
        return piHap;
    }

        /**
    * Find initial values of parameters representing piHap frequencies of haplotypes
    */
    private void initialiseHapParams(File fileNameFile, boolean initialFreqSet, List<String> initialFreqVector) {
        double remaining = 1.0;
        initialHapParams = new double[nTimePoints][nHaplo-1];  // Parameters encoding haplotype frequencies
        for (int iHaplo = 0; iHaplo < nHaplo-1; iHaplo++) {
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                initialHapParams[iTimePoint][iHaplo] = (1.0 + 0.1 * (random.nextDouble() - 0.5)) /(nHaplo-iHaplo);
            }
        }
        
        nActive = new int[nTimePoints];
        Arrays.fill(nActive, nHaplo-1);
        mapping = new int[nTimePoints][nHaplo-1];                       // Mapping of parameters to deal with fixed frequencies
        for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
            for (int iHaplo = 0; iHaplo < nHaplo-1; iHaplo++) {
                mapping[iTimePoint][iHaplo] = iHaplo;
            }
        }
        fixedValues = new double[nTimePoints][nHaplo-1];


        for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {  // read in data files
            System.out.println("initialFreqSet = " + initialFreqSet);
            if (initialFreqSet) {
                String[] words = initialFreqVector.get(iTimePoint).split("\\t");
                System.out.println(words.length + "\t" + nHaplo);
                if (words.length > nHaplo) {
                    double[] initialHapFreq = new double[nHaplo];
                    Arrays.fill(mapping[iTimePoint], -1);
                    int iPoint = -1;
                    boolean fixed = false;
                    for (int iHaplo = 0; iHaplo < nHaplo-1; iHaplo++) {
                        initialHapFreq[iHaplo] = Double.parseDouble(words[iHaplo+1]);
                        if (initialHapFreq[iHaplo] >= 1.) {
                            initialHapFreq[iHaplo] = 1.0;
                            Arrays.fill(fixedValues[iTimePoint], 0.0);
                            Arrays.fill(mapping[iTimePoint], -1);
                            fixedValues[iTimePoint][iHaplo] = 1.0;
                            iPoint = -1;
                            fixed = true;
                            break;
                        } else if (initialHapFreq[iHaplo] <= 0.) {
                            initialHapFreq[iHaplo] = 0.0;
                            fixedValues[iTimePoint][iHaplo] = 0.0;
                        } else {
                            iPoint++;
                            mapping[iTimePoint][iPoint] = iHaplo;
                        }
                    }
                    initialHapFreq[nHaplo-1] = Double.parseDouble(words[nHaplo]);
                    if (initialHapFreq[nHaplo-1] >= 1.) {
                        initialHapFreq[nHaplo-1] = 1.0;
                        Arrays.fill(fixedValues[iTimePoint], 0.0); 
                        Arrays.fill(mapping[iTimePoint], -1);
                        iPoint = -1;
                        fixed = true;
                    } else if (!fixed && iPoint > -1 && initialHapFreq[nHaplo-1] <= 0.) {
                        initialHapFreq[nHaplo-1] = 0.0;
                        System.out.println("iPoint = " + iPoint);
                        System.out.println("mapping[iTimePoint][iPoint] = " + mapping[iTimePoint][iPoint]);
                        fixedValues[iTimePoint][mapping[iTimePoint][iPoint]] = 1.0;
                        mapping[iTimePoint][iPoint] = -1;
                        iPoint--;
                    }
                    nActive[iTimePoint] = iPoint+1;
                    if (true) {
                        System.out.println(initialFreqVector.get(iTimePoint));
                        System.out.println("nActive: " + nActive[iTimePoint] + 
                                "\tfixed values: " + Arrays.toString(fixedValues[iTimePoint])
                                + "\tmappings: " + Arrays.toString(Arrays.copyOf(mapping[iTimePoint], nActive[iTimePoint])));
                    }
                    double[] scaledFreq = new double[nHaplo-1];
                    double remainder = 1.0;
                    for (int iHaplo = 0; iHaplo < nHaplo - 1; iHaplo++) {
                        scaledFreq[iHaplo] = initialHapFreq[iHaplo]/remainder;
                        remainder -= initialHapFreq[iHaplo];
                    }
                    Arrays.fill(initialHapParams[iTimePoint], 0.0);
                    for (int iParam = 0; iParam < nActive[iTimePoint]; iParam++) {
                        initialHapParams[iTimePoint][iParam] = scaledFreq[mapping[iTimePoint][iParam]];
                    }
//                    System.out.println(Arrays.toString(scaledFreq));
                }
            }
//            System.out.println(iTimePoint + "\t" + Arrays.toString(Arrays.copyOf(initialHapParams[iTimePoint], nActive[iTimePoint])));
        }     
    }
    
        
    double[][] getInitialHapParams() {
        return initialHapParams;
    }
    
    int getNTimePoints() {
        return nTimePoints;
    }


    public int getSiteCount() {
        return siteCount;
    }
    
    public int[] getNActive() {
        return nActive;
    }
    
    
}
    


