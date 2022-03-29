/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package cluster_rg;

import java.util.Arrays;

/**
 *
 * @author rgoldst
 */
public class Assignment {
    
    int[] assign = null;
    boolean[] presentBase = new boolean[4];
    int nPresent = 0;

    private int nHaplo = 0;
    private int nAbsent = 0;
    private int nTimePoints = 0;

    private final GammaCalc gammaCalc;
    private final boolean verbose;
    
    Assignment(int iAssign, int nHaplo, GammaCalc gammaCalc, boolean verbose) {
        this.gammaCalc = gammaCalc;
        this.nHaplo = nHaplo;
        this.verbose = verbose;
        assign = new int[nHaplo];
        for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {    // Loop over possible haplotypes
            assign[iHaplo] = (iAssign / pow(Constants.MAX_BASES, iHaplo)) % (Constants.MAX_BASES);
            presentBase[assign[iHaplo]] = true;
        }       
        for (int iBase = 0; iBase < 4; iBase++) {
            if (presentBase[iBase]) {
                nPresent++;
            }
        }
        nAbsent = 4 - nPresent;
        if (this.verbose) {
            System.out.println(iAssign + "\t" + Arrays.toString(assign));
        }
    }

    double computeAssignmentLogLikelihood(int iTimePoint, int[][] strandReads, int[] totStrand,
            boolean siteConserved, double alpha0, double alphaE, double[] piHap) {
        double[] piNuc = new double[4];
        double[] alphaObs = new double[4];
        double sumAlphaObs = 0.0;
        for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
            piNuc[assign[iHaplo]] += piHap[iHaplo];
        }
        for (int iBase = 0; iBase < 4; iBase++) {
            alphaObs[iBase] = piNuc[iBase] * alpha0 +  alphaE;
            sumAlphaObs += alphaObs[iBase];
        }
        double[] logLikelihoodStrand = new double[2];
        double g1 = this.gammaCalc.logGamma(sumAlphaObs);
        for (int iStrand = 0; iStrand < 2; iStrand++) {
            logLikelihoodStrand[iStrand] = g1 - this.gammaCalc.logGamma(sumAlphaObs + totStrand[iStrand]);
            for (int iBase = 0; iBase < 4; iBase++) {
                if (strandReads[iStrand][iBase] > 0) {
                    logLikelihoodStrand[iStrand] += this.gammaCalc.logGamma(alphaObs[iBase] + strandReads[iStrand][iBase])
                            - this.gammaCalc.logGamma(alphaObs[iBase]);
                }
            }          
        }
        double logFitness = logLikelihoodStrand[0] + logLikelihoodStrand[1];
        return logFitness;
    }
        
//    double computeExpectedLikelihood(int iTimePoint, int[][] strandReads, int[] reads, int[] totStrand, boolean siteConserved , boolean printStuff) {
//        double[][] logLikelihoodStrand = new double[2][2];
//        int[][] sampleStrandReads = new int[2][4];
//        int[] sampleReads = new int[4];
//        int[] sampleTotStrand = new int[2];
//        double[] baseProb = new double[4];
//        for (int iStrand = 0; iStrand < 2; iStrand++) {
//            for (int iBase = 0; iBase < 4; iBase++) {
//                baseProb[iBase] = 
//                        (currentPiNuc[iTimePoint][iBase] * currentAlpha0 + (1.0 - currentPiNuc[iTimePoint][iBase]) * currentAlphaE)
//                        / (currentAlpha0 + 3.0 * currentAlphaE);
//                double stdDev = Math.sqrt(totStrand[iStrand] * baseProb[iBase] * (1.0 - baseProb[iBase]));
//                double expBaseCount = totStrand[iStrand] * baseProb[iBase];
//                sampleStrandReads[iStrand][iBase] = Math.round(Math.round(
//                        expBaseCount + Math.pow(-1.0, iBase+iStrand) * stdDev));
//                sampleStrandReads[iStrand][iBase] 
//                        = Math.max(0, Math.min(totStrand[iStrand], sampleStrandReads[iStrand][iBase]));
//                sampleReads[iBase] += sampleStrandReads[iStrand][iBase];
//                sampleTotStrand[iStrand] += sampleStrandReads[iStrand][iBase]; 
//            }
//        }
//        return computeAssignmentLogLikelihood(iTimePoint, sampleStrandReads, sampleReads, sampleTotStrand, false);
//    }
    
    int pow (int a, int b) {  // Computes powers
        if ( b == 0)     return 1;
        if ( b == 1)     return a;
        if (b%2 == 0)    return     pow ( a * a, b/2); //even a=(a^2)^b/2
        else             return a * pow ( a * a, (b-1)/2); //odd  a=a*(a^2)^b/2
    }

    
}
