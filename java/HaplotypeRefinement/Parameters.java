import java.util.ArrayList;
import java.util.Arrays;

/**
 * Contains set of parameters of a model and methods to alter the 
 * form of that model through recombining, expanding, and contracting
 * 
 * @author rgoldst
 */
public class Parameters {

    double totalLogLikelihood = 0.;         // Log likelihood of counts
    double totPenalty = 0.;                 // Penalty for parameters
    double[][][] probBase = null;           // Base frequencies for each haplotype
    double[] hapFreq = null;                // Haplotype frequencies
    int[][] iMaxBase = null;                // Most common base
    int nHaplo = 0;                         // Number of haplotypes
    double adjLogLikelihood = 0.;           // Log likelihood corrected for number of parameters
    ArrayList<Integer>[] polySites;         // List of variable sites
    int nActive = 0;                        // Number of active haplotypes with non-negligible frequencies

    /** 
     * 
     * @param totalLogLikelihood
     * @param totPenalty
     * @param hapFreq
     * @param probBase
     * @param iMaxBase
     * @param polySites 
     */
    Parameters(double totalLogLikelihood, double totPenalty, 
            double[] hapFreq, double[][][] probBase,
            int[][] iMaxBase, ArrayList<Integer>[] polySites) {
        this.totalLogLikelihood = totalLogLikelihood;
        this.totPenalty = totPenalty;
        this.adjLogLikelihood = totalLogLikelihood - totPenalty;
        this.hapFreq = hapFreq;
        this.probBase = probBase;
        this.iMaxBase = iMaxBase;
        this.nHaplo = hapFreq.length;
        this.polySites = polySites;
        nActive = 0;
        for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
            if (hapFreq[iHaplo] > 1.0E-10) {
                nActive++;
            }
        }
    } 
      
    /**
     * Make a new model by recombining haplotypes, and return new Parameters file
     * 
     * @param nSites
     * @param activeSites
     * @return 
     */
    Parameters recombineHaplo(int nSites, ArrayList<Integer> activeSites) {
        double[] newHapFreq = Arrays.copyOf(hapFreq, nHaplo);                               // Transfer current state: haplotype frequencies
        double[][][] newProbBase = new double[nHaplo][nSites][4];
        for (int iSite : activeSites) {                                                     // Transfer current state: base probabilities
            for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                newProbBase[iHaplo][iSite] = Arrays.copyOf(probBase[iHaplo][iSite],4);
            }
        }

        boolean ok = false;                                                                 // Boolean: did recombination step work?
        int iTry = 0;
        int iHaplo = 0;                                                                     // Two haplotypes to recombine
        int jHaplo = 0;
        int lBreak = 0;
        int rBreak = 0;
        int ranForm = 0;
        while (!ok) {
            iHaplo = RefineHaplotypes.random.nextInt(nHaplo);                               // Choose two haplotypes at random
            jHaplo = (iHaplo + RefineHaplotypes.random.nextInt(nHaplo-1) + 1)%nHaplo;
            ranForm = RefineHaplotypes.random.nextInt(4);                            // Choose between recombination switch or gene conversion
            
            int iBreak = RefineHaplotypes.random.nextInt(nSites);                           // Set left and right breaks with length selected by normal distribution
            int jBreak = iBreak + Math.round(Math.round(RefineHaplotypes.nd.sample()));
            lBreak = Math.min(iBreak, jBreak);
            rBreak = Math.max(iBreak, jBreak);

            int nSwitch = 0;
            if (hapFreq[iHaplo] > 0.01 && hapFreq[jHaplo] > 0.01)  {                        // Don't bother with very low frequency haplotyes
                for (int iSite : activeSites) {
                    if (iSite > lBreak && iSite < rBreak) {                                 // Within break region
                        if (ranForm == 2) {                                                 // Overwrite i with j in this region
                            newProbBase[iHaplo][iSite] = Arrays.copyOf(probBase[jHaplo][iSite],4);
                            newProbBase[jHaplo][iSite] = Arrays.copyOf(probBase[jHaplo][iSite],4);
                        } else if (ranForm == 3) {                                          // Overwrite j with i in this region
                            newProbBase[iHaplo][iSite] = Arrays.copyOf(probBase[iHaplo][iSite],4);
                            newProbBase[jHaplo][iSite] = Arrays.copyOf(probBase[iHaplo][iSite],4);
                        } else {                                                            // Switch i and j in this region
                            newProbBase[iHaplo][iSite] = Arrays.copyOf(probBase[jHaplo][iSite],4);
                            newProbBase[jHaplo][iSite] = Arrays.copyOf(probBase[iHaplo][iSite],4);                            
                        }
                        nSwitch++;
                    }
                }
            }
            iTry++;
            ok = (nSwitch > 0 || iTry > 10);                                                // If non-conserved sites actually switched or tried too many times declare success
        }
        if (ranForm == 2) {
            System.out.println("Overwriting " + iHaplo + " with " + jHaplo + " in region " + lBreak + "-" + rBreak);
            RefineHaplotypes.logWriter.println("Overwriting " + iHaplo + " with " + jHaplo + " in region " + lBreak + "-" + rBreak);
        } else if (ranForm == 3) {
            System.out.println("Overwriting " + jHaplo + " with " + iHaplo + " in region " + lBreak + "-" + rBreak);
            RefineHaplotypes.logWriter.println("Overwriting " + jHaplo + " with " + iHaplo + " in region " + lBreak + "-" + rBreak);

        } else {
            System.out.println("Recombining haplotypes " + iHaplo + " and " + jHaplo + " in region " + lBreak + "-" + rBreak);
            RefineHaplotypes.logWriter.println("Recombining haplotypes " + iHaplo + " and " + jHaplo + " in region " + lBreak + "-" + rBreak);
        }
        Parameters newParameters = new Parameters(0.0, 0.0, newHapFreq, newProbBase, 
            null, null);
        return newParameters;
            
    }
   
    /**
     * Create a new model by compressing two haplotypes into one, and return new Parameters file
     * 
     * @param iHaplo
     * @param jHaplo
     * @param nSites
     * @return 
     */
    Parameters compressHaplo(int iHaplo, int jHaplo, int nSites) {
        System.out.println("Combining haplotypes " + iHaplo + " and " + jHaplo);
        RefineHaplotypes.logWriter.println("Combining haplotypes " + iHaplo + " and " + jHaplo);      
        double[] contractedHapFreq = new double[nHaplo-1];                                              // Generate new Parameters file with one fewer haplotype
        double[][][] contractedProbBase = new double[nHaplo-1][nSites][4];                              // Reduced base frequencies
        int iPoint = 0;                                                                                 // Copy unaffected parameters into first n-2 = n'-1 haplotypes
        for (int kHaplo = 0; kHaplo < nHaplo; kHaplo++) {
            if (kHaplo != iHaplo && kHaplo != jHaplo) {
                contractedHapFreq[iPoint] = hapFreq[kHaplo];
                for (int iSite = 0; iSite < nSites; iSite++) {
                    contractedProbBase[iPoint][iSite] = Arrays.copyOf(probBase[kHaplo][iSite], 4);
                }  
                iPoint++;
            }
        }

        contractedHapFreq[nHaplo-2] = hapFreq[iHaplo] + hapFreq[jHaplo];                                // New haplotype has frequency = sum of parent haplotypes
        for (int iSite = 0; iSite < nSites; iSite++) {                                                  // And base frequencies are averages of two parents
            for (int iBase = 0; iBase < 4; iBase++) {
                contractedProbBase[nHaplo-2][iSite][iBase] = 
                        0.5 * (probBase[iHaplo][iSite][iBase] + probBase[jHaplo][iSite][iBase]);
            }
        }      
        int contractedNHaplo = nHaplo-1;
            System.out.println("Combining haplotypes " + iHaplo + " and " + jHaplo);
            RefineHaplotypes.logWriter.println("Combining haplotypes " + iHaplo + " and " + jHaplo);
            
        Parameters newParameters = new Parameters(0.0, 0.0, contractedHapFreq, contractedProbBase, 
            null, null);
        return newParameters;      
    }
    
    /**
     * Create a new model by expanding a haplotype, and return new Parameters file
     * 
     * @param expandHaplo
     * @param nSites
     * @return 
     */
    Parameters expandHaplo(int expandHaplo, int nSites) {           
        double[] expandedHapFreq = new double[nHaplo+1];                                                    // Prepare new set of parameters with one more haplotype
        double[][][] expandedProbBase = new double[nHaplo+1][nSites][4];

        // Copy parameters for expanded set
        for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {                                                   // Copy current ones into next
            expandedHapFreq[iHaplo] = hapFreq[iHaplo];
            for (int iSite = 0; iSite < nSites; iSite++) {
                expandedProbBase[iHaplo][iSite] = Arrays.copyOf(probBase[iHaplo][iSite], 4);
            }
        }
        for (int iSite = 0; iSite < nSites; iSite++) {
            expandedProbBase[nHaplo][iSite] = Arrays.copyOf(probBase[expandHaplo][iSite], 4);
        }  
        double avgFrac = 0.0;

        // Provide initial base probabilities for expanded pairs
        for (int iSite : polySites[expandHaplo]){
            int maxBase = -1;
            double maxBaseProb = -1.0;
            for (int iBase = 0; iBase < 4; iBase++) {
                if (probBase[expandHaplo][iSite][iBase] > maxBaseProb) {
                    maxBaseProb = probBase[expandHaplo][iSite][iBase];
                    maxBase = iBase;
                }
                expandedProbBase[expandHaplo][iSite][iBase] = 0.2 * probBase[expandHaplo][iSite][iBase];    // New haplo gets most frequent base of old
                expandedProbBase[nHaplo][iSite][iBase] = 0.8 * probBase[expandHaplo][iSite][iBase];         // Revised old gets expanded others
            }
            expandedProbBase[expandHaplo][iSite][maxBase] = 0.8 * probBase[expandHaplo][iSite][maxBase];
            expandedProbBase[nHaplo][iSite][maxBase] = 0.2 * probBase[expandHaplo][iSite][maxBase];

            // Make sure it is adequately normalised
            double summ = expandedProbBase[expandHaplo][iSite][0]+expandedProbBase[expandHaplo][iSite][1]
                +expandedProbBase[expandHaplo][iSite][2]+expandedProbBase[expandHaplo][iSite][3];
            expandedProbBase[expandHaplo][iSite][0]/=summ;
            expandedProbBase[expandHaplo][iSite][1]/=summ;
            expandedProbBase[expandHaplo][iSite][2]/=summ;
            expandedProbBase[expandHaplo][iSite][3]/=summ;
            summ = expandedProbBase[nHaplo][iSite][0]+expandedProbBase[nHaplo][iSite][1]
                +expandedProbBase[nHaplo][iSite][2]+expandedProbBase[nHaplo][iSite][3];
            expandedProbBase[nHaplo][iSite][0]/=summ;
            expandedProbBase[nHaplo][iSite][1]/=summ;
            expandedProbBase[nHaplo][iSite][2]/=summ;
            expandedProbBase[nHaplo][iSite][3]/=summ;
            avgFrac += maxBaseProb;
        }            

        avgFrac /= polySites[expandHaplo].size();
        expandedHapFreq[expandHaplo] = avgFrac*hapFreq[expandHaplo];
        expandedHapFreq[nHaplo] = (1.0-avgFrac)*hapFreq[expandHaplo];
        int expandedNHaplo = nHaplo+1;
        System.out.println("\nExpanding haplotype " + expandHaplo);
        RefineHaplotypes.logWriter.println("\nExpanding haplotype " + expandHaplo);     
        Parameters newParameters = new Parameters(0.0, 0.0, expandedHapFreq, expandedProbBase, 
            null, null);

        return newParameters; 
    }
    


    
}
