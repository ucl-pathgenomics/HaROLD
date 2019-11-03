/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package refineHaplotypes;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 *
 * @author rgoldst
 */
public class Parameters {

    double totalLogLikelihood = 0.;
    double totPenalty = 0.;
    double[][][] probBase = null;
    double[] hapFreq = null;
    int[][] iMaxBase = null;
    int nHaplo = 0;
    double adjLogLikelihood = 0.;
    ArrayList<Integer>[] polySites;
    int nActive = 0;

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
        
    Parameters recombineHaplo(int nSites, ArrayList<Integer> activeSites) {
        double[] newHapFreq = Arrays.copyOf(hapFreq, nHaplo);
        double[][][] newProbBase = new double[nHaplo][nSites][4];
        for (int iSite : activeSites) {
            for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                newProbBase[iHaplo][iSite] = Arrays.copyOf(probBase[iHaplo][iSite],4);
            }
        }

        boolean ok = false;
        int iTry = 0;
        int iHaplo = 0;
        int jHaplo = 0;
        while (!ok) {
            iHaplo = RefineHaplotypes.random.nextInt(nHaplo);
            jHaplo = (iHaplo + RefineHaplotypes.random.nextInt(nHaplo-1) + 1)%nHaplo;
            double ranForm = RefineHaplotypes.random.nextInt(4);
            
            int iBreak = RefineHaplotypes.random.nextInt(nSites);
            int jBreak = iBreak + Math.round(Math.round(RefineHaplotypes.nd.sample()));
            int lBreak = Math.min(iBreak, jBreak);
            int rBreak = Math.max(iBreak, jBreak);

            int nSwitch = 0;
            if (hapFreq[iHaplo] > 0.01 && hapFreq[jHaplo] > 0.01)  {
                for (int iSite : activeSites) {
                    if (iSite > lBreak && iSite < rBreak) {
                        if (ranForm == 2) {
                            newProbBase[iHaplo][iSite] = Arrays.copyOf(probBase[jHaplo][iSite],4);
                            newProbBase[jHaplo][iSite] = Arrays.copyOf(probBase[jHaplo][iSite],4);
                        } else if (ranForm == 3) {
                            newProbBase[iHaplo][iSite] = Arrays.copyOf(probBase[iHaplo][iSite],4);
                            newProbBase[jHaplo][iSite] = Arrays.copyOf(probBase[iHaplo][iSite],4);
                        } else {
                            newProbBase[iHaplo][iSite] = Arrays.copyOf(probBase[jHaplo][iSite],4);
                            newProbBase[jHaplo][iSite] = Arrays.copyOf(probBase[iHaplo][iSite],4);                            
                        }
                        nSwitch++;
                    }
                }
            }
            iTry++;
            ok = (nSwitch > 0 || iTry > 10);
        }


        System.out.println("Recombining haplotypes " + iHaplo + " and " + jHaplo);
        RefineHaplotypes.logWriter.println("Combining haplotypes " + iHaplo + " and " + jHaplo);      
        Parameters newParameters = new Parameters(0.0, 0.0, newHapFreq, newProbBase, 
            null, null);

        return newParameters;                       // Did not expand
            
    }

    Parameters compressHaplo(int iHaplo, int jHaplo, int nSites) {
        System.out.println("Combining haplotypes " + iHaplo + " and " + jHaplo);
        RefineHaplotypes.logWriter.println("Combining haplotypes " + iHaplo + " and " + jHaplo);      
        double[] contractedHapFreq = new double[nHaplo-1];
        double[][][] contractedProbBase = new double[nHaplo-1][nSites][4];
        int iPoint = 0;
        for (int kHaplo = 0; kHaplo < nHaplo; kHaplo++) {
            if (kHaplo != iHaplo && kHaplo != jHaplo) {
                contractedHapFreq[iPoint] = hapFreq[kHaplo];
                for (int iSite = 0; iSite < nSites; iSite++) {
                    contractedProbBase[iPoint][iSite] = Arrays.copyOf(probBase[kHaplo][iSite], 4);
                }  
                iPoint++;
            }
        }

        contractedHapFreq[nHaplo-2] = hapFreq[iHaplo] + hapFreq[jHaplo];
        for (int iSite = 0; iSite < nSites; iSite++) {
            for (int iBase = 0; iBase < 4; iBase++) {
                contractedProbBase[nHaplo-2][iSite][iBase] = 
                        0.5 * (probBase[iHaplo][iSite][iBase] + probBase[jHaplo][iSite][iBase]);
            }
        }      
        int contractedNHaplo = nHaplo-1;
        Parameters newParameters = new Parameters(0.0, 0.0, contractedHapFreq, contractedProbBase, 
            null, null);
        return newParameters;      
    }
    
    
    Parameters expandHaplo(int expandHaplo, int nSites) {
        System.out.println("\nExpanding haplotype " + expandHaplo);
        RefineHaplotypes.logWriter.println("\nExpanding haplotype " + expandHaplo);            
        double[] expandedHapFreq = new double[nHaplo+1];
        double[][][] expandedProbBase = new double[nHaplo+1][nSites][4];

        // Copy parameters for expanded set
        for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
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
                expandedProbBase[expandHaplo][iSite][iBase] = 0.2 * probBase[expandHaplo][iSite][iBase];
                expandedProbBase[nHaplo][iSite][iBase] = 0.8 * probBase[expandHaplo][iSite][iBase];
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

        Parameters newParameters = new Parameters(0.0, 0.0, expandedHapFreq, expandedProbBase, 
            null, null);

        return newParameters; 
    }
    


    
}
