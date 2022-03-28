/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package cluster_rg;

import java.util.ArrayList;
import java.util.Arrays;


/**
 * Container for the reads at a single site
 * 
 * @author rgoldst
 */
public class Site {   
    int iSite;
    double[] estProbDiffBases = new double[5];
    int conservedBase = -9;
    boolean siteConserved = false;

    private int nTimePoints;
    private int nHaplo = 0;
    private int nBases = 0;
    private ArrayList<Assignment> assignmentVector;
    private ArrayList<Assignment> localAssignmentVector = new ArrayList<>();
    private double[] probAssignment = null;
    private double[] priorProb = new double[2];
    private int nAssignments = 0;
    private int[] nActive;
    private String[] baseString = {"A", "C", "G", "T"};
    private int bestGuessNPresent = -1;
    double expectedLikelihood;
    double totalLogLikelihood;

     int[][][] strandReads = null; // [tp][strand][base] top two sets of reads on each strand
    private int[][] totStrand = null; // [tp][strand] number of reads on each strand
    private int[][] reads = null; // [tp][base]
    private int[] totReads = null; // [tp]
    private int[] timePointConservedBase = null;

    private boolean siteActive = false;
    private int nPresentBase = 0;   // number of present bases

    private boolean[] presentBase = new boolean[4];
    private boolean[] timePointConserved = null;
    private boolean[] timePointHasData = null;

    private final GammaCalc gammaCalc;

    Site(int iSite, int nTimePoints, int[] nActive, int nHaplo, ArrayList<Assignment> assignmentVector, GammaCalc gammaCalc) {
        this.gammaCalc = gammaCalc;
        this.iSite = iSite;
        this.nTimePoints = nTimePoints;
        this.nActive = nActive;
        this.nHaplo = nHaplo;
        this.assignmentVector = assignmentVector;
        nAssignments = assignmentVector.size();
        strandReads = new int[nTimePoints][2][4]; // [tp][strand][base] top two sets of reads on each strand
        totStrand = new int[nTimePoints][2]; // [tp][strand] number of reads on each strand
        reads = new int[nTimePoints][4]; // [tp][base]
        totReads = new int[nTimePoints]; // [tp]
        timePointConserved = new boolean[nTimePoints];
        timePointHasData = new boolean[nTimePoints];
        timePointConservedBase = new int[nTimePoints];
        for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
            timePointHasData[iTimePoint] = false; // assume timepoint does not have data
            timePointConserved[iTimePoint] = false; // assume timepoint is not conserved
        }
    }

    // Compute probAssignment[], representing probability of a given assignment of bases to haplotypes at that site
    // Return likelihood
    double assignHaplotypes(double alpha0, double alphaE, double lGA0, double lGA03AE, 
            double[] priors, double[][] piHap) {
        double logLikelihood = 0.0;
        if (siteConserved) {
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                for (int iStrand = 0; iStrand < 2; iStrand++) {
                    logLikelihood += lGA03AE - lGA0
                            - this.gammaCalc.logGamma(alpha0 + 3.0 * alphaE + totStrand[iTimePoint][iStrand])
                            + this.gammaCalc.logGamma(alpha0 + totStrand[iTimePoint][iStrand]);
                }
            }
            bestGuessNPresent = 1;
            return (priors[1] + logLikelihood);
        }
            
        probAssignment = new double[localAssignmentVector.size()];
        double[] logLikelihoodAssign = new double[localAssignmentVector.size()];
        double sumProb = 0.0;
        int bestAssign = -999;
        double bestAssignVal = -1.0E20;
        for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
            logLikelihoodAssign[iAssign] = priors[localAssignmentVector.get(iAssign).nPresent];
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                logLikelihoodAssign[iAssign] += 
                        + localAssignmentVector.get(iAssign).computeAssignmentLogLikelihood(iTimePoint, strandReads[iTimePoint], 
                                totStrand[iTimePoint], siteConserved, alpha0, alphaE, piHap[iTimePoint]);
            }
            if (logLikelihoodAssign[iAssign] > bestAssignVal) {
                bestAssignVal = logLikelihoodAssign[iAssign];
                bestAssign = iAssign;
            }
        }
        if (bestAssign < 0) {
            System.out.println("LogLiklihoodAssign: " + Arrays.toString(logLikelihoodAssign));
            for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
                logLikelihoodAssign[iAssign] = priors[localAssignmentVector.get(iAssign).nPresent];
                System.out.println(iAssign + "\t" + logLikelihoodAssign[iAssign]);
                for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                    logLikelihoodAssign[iAssign] += 
                            + localAssignmentVector.get(iAssign).computeAssignmentLogLikelihood(iTimePoint, strandReads[iTimePoint], 
                                    totStrand[iTimePoint], siteConserved, alpha0, alphaE, piHap[iTimePoint]);
                    System.out.println(iTimePoint + "\t" + Arrays.toString(strandReads[iTimePoint][0]) + Arrays.toString(strandReads[iTimePoint][1]) 
                            + "\t" + Arrays.toString(totStrand[iTimePoint])
                            + "\t" + alpha0 + "\t" + alphaE + "\t" + Arrays.toString(piHap[iTimePoint]));
                    System.out.println("ComputeAssignLogLik\t" + localAssignmentVector.get(iAssign).computeAssignmentLogLikelihood(iTimePoint, strandReads[iTimePoint], 
                                    totStrand[iTimePoint], siteConserved, alpha0, alphaE, piHap[iTimePoint]));
                }
            }
            System.exit(1);
        }
        bestGuessNPresent = localAssignmentVector.get(bestAssign).nPresent;
        
        for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
            probAssignment[iAssign] = Math.exp(logLikelihoodAssign[iAssign]-bestAssignVal);
            sumProb += probAssignment[iAssign];
            logLikelihood += probAssignment[iAssign];
        }
        for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
            probAssignment[iAssign] /= sumProb;
        }
        totalLogLikelihood = bestAssignVal + Math.log(logLikelihood);
        return totalLogLikelihood;
    }
       
    double computeSiteLogLikelihood(double alpha0, double alphaE, 
            double lGA0, double lGA03AE, double[] priors, double[][] piHap) {
        double logLikelihood = 0.0; 
        if (siteConserved) {
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                for (int iStrand = 0; iStrand < 2; iStrand++) {
                    logLikelihood += lGA03AE - lGA0
                            - this.gammaCalc.logGamma(alpha0 + 3.0 * alphaE + totStrand[iTimePoint][iStrand])
                            + this.gammaCalc.logGamma(alpha0 + totStrand[iTimePoint][iStrand]) ;
                }
            }
            double totalLogLikelihood = priors[1] + logLikelihood;
            return totalLogLikelihood;
        }   
        // As the sum is over all timepoints for each assignment
        double bestAssignVal = -1.0E10;
        double[] logLikelihoodAssign = new double[localAssignmentVector.size()];
        for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
            if (probAssignment[iAssign] > 0.001) {
                logLikelihoodAssign[iAssign] = priors[localAssignmentVector.get(iAssign).nPresent];
                for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                    logLikelihoodAssign[iAssign] += 
                            localAssignmentVector.get(iAssign).computeAssignmentLogLikelihood(iTimePoint, strandReads[iTimePoint], 
                                    totStrand[iTimePoint], siteConserved, alpha0, alphaE, piHap[iTimePoint]);
                }
                if (Double.isInfinite(logLikelihoodAssign[iAssign]) || logLikelihoodAssign[iAssign] < -1000.0) {
                    logLikelihoodAssign[iAssign] = -1000.0;
                }
                bestAssignVal = Math.max(bestAssignVal, logLikelihoodAssign[iAssign]);
            }
        }
        for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
            if (probAssignment[iAssign] > 0.001) {
                logLikelihood += Math.exp(logLikelihoodAssign[iAssign]-bestAssignVal);
            }
        }
        totalLogLikelihood = bestAssignVal + Math.log(logLikelihood);
        
        return totalLogLikelihood;
    }
    

    // Compute log likelihood for a given time point, using already calculated probability of assignment
    double computeSiteTimePointLogLikelihood(int iTimePoint, double alpha0, 
            double alphaE, double lGA0, double lGA03AE, double[] priors, double[] piHap) {
        double logLikelihood = 0.0;
        if (siteConserved) {
            for (int iStrand = 0; iStrand < 2; iStrand++) {
                logLikelihood += lGA03AE - lGA0
                        - this.gammaCalc.logGamma(alpha0 + 3.0 * alphaE + totStrand[iTimePoint][iStrand])
                        + this.gammaCalc.logGamma(alpha0 + totStrand[iTimePoint][iStrand]);
            }
            return logLikelihood;
        }
        double[] logLikelihoodAssign = new double[localAssignmentVector.size()];
        double bestAssignVal = -1.0E20;
        // Note because we are dealing with a single point we include probAssignment but not priors
        for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
            if (probAssignment[iAssign] > 0.001) {
                logLikelihoodAssign[iAssign] =
                            localAssignmentVector.get(iAssign).computeAssignmentLogLikelihood(iTimePoint, strandReads[iTimePoint], 
                                    totStrand[iTimePoint], siteConserved, alpha0, alphaE, piHap);
                bestAssignVal = Math.max(bestAssignVal, logLikelihoodAssign[iAssign]);
            }
        }
        for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
            if (probAssignment[iAssign] > 0.001) {
                logLikelihood += probAssignment[iAssign] * Math.exp(logLikelihoodAssign[iAssign]-bestAssignVal);
            }
        }
        double siteLogLikelihood = bestAssignVal + Math.log(logLikelihood);
        return siteLogLikelihood;
    }
    
//    void computeExpectedLikelihood(double alpha0, double alphaE, double lGA0, double lGA03AE, double[] priors, ) {
//        double logLikelihood = 0.0;
//        double[] logLikelihoodAssign = new double[localAssignmentVector.size()];
//        double bestAssignVal = -1.0E20;
//        for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
//            logLikelihoodAssign[iAssign] = priors[localAssignmentVector.get(iAssign).nPresent];
//            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
//                logLikelihoodAssign[iAssign] += 
//                        + localAssignmentVector.get(iAssign).computeExpectedLikelihood(iTimePoint, strandReads[iTimePoint], 
//                                reads[iTimePoint], totStrand[iTimePoint], siteConserved, false);
//            }
//            if (logLikelihoodAssign[iAssign] > bestAssignVal) {
//                bestAssignVal = logLikelihoodAssign[iAssign];
//            }
//        }
//        for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
//            logLikelihood += Math.exp(logLikelihoodAssign[iAssign]-bestAssignVal);
//        }
//        expectedLikelihood = bestAssignVal + Math.log(logLikelihood);
//    }

    
    double[][] getProbBase() {
        double[][] expectedFreq = new double[nHaplo][4];
        if (siteConserved) {
            for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                expectedFreq[iHaplo][conservedBase] = 1.0;
            }
        } else {
            for (int iAssign = 0; iAssign < localAssignmentVector.size(); iAssign++) {
                Assignment assignment = localAssignmentVector.get(iAssign);
                for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                    expectedFreq[iHaplo][assignment.assign[iHaplo]] += probAssignment[iAssign];
                }
            }
        }
        if (!siteConserved) { 
            if (smellTest(false)) {
                for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                    for (int iBase = 0; iBase < 4; iBase++) {
                        expectedFreq[iHaplo][iBase] = 0.0;
                    }
                }
            }
        }
        return expectedFreq;
    }

      
    void addTimePoint(int iTimePoint, String line) { 
        String[] words = line.split(",");
        for (int iBase = 0; iBase < 4; iBase++) {   // compute various sums of reads
            
            if (words.length > 10) {
                for (int iStrand = 0; iStrand < 2; iStrand++) {
                    strandReads[iTimePoint][iStrand][iBase]=Integer.parseInt(words[(2*iBase)+iStrand+3]);
                    reads[iTimePoint][iBase] += strandReads[iTimePoint][iStrand][iBase];
                }
                for (int iStrand = 0; iStrand < 2; iStrand++) {
                    totStrand[iTimePoint][iStrand] += strandReads[iTimePoint][iStrand][iBase];
                    totReads[iTimePoint] += strandReads[iTimePoint][iStrand][iBase];
                }
            } else {
                strandReads[iTimePoint][0][iBase]=Integer.parseInt(words[(iBase)+3]);
                reads[iTimePoint][iBase] += strandReads[iTimePoint][0][iBase];
                totStrand[iTimePoint][0] += strandReads[iTimePoint][0][iBase];
                totReads[iTimePoint] += strandReads[iTimePoint][0][iBase];          
            }
            
            if (reads[iTimePoint][iBase] > 0) {
                conservedBase = iBase;
                presentBase[iBase] = true;
            }
        }
        if (totReads[iTimePoint] > 0) {   // has data
            timePointHasData[iTimePoint] = true;
        }
        if (Math.max(Math.max(reads[iTimePoint][0],reads[iTimePoint][1]),
                Math.max(reads[iTimePoint][2],reads[iTimePoint][3])) == totReads[iTimePoint]) {  // timePointConserved site
            timePointConserved[iTimePoint] = true;
            for (int iBase = 0; iBase < 4; iBase++) {
                if (reads[iTimePoint][iBase] == totReads[iTimePoint]) {
                    timePointConservedBase[iTimePoint] = iBase;
                }
            }
        }
    }

    
    boolean isActive() {
        siteActive = false;
        for (int iBase = 0; iBase < 4; iBase++) {
            if (presentBase[iBase]) {
                nPresentBase++;
                siteActive = true;
            }
        }
        siteConserved = (nPresentBase == 1);
        for (Assignment assignment : assignmentVector) {
            boolean addThis = true;
            for (int iBase = 0; iBase < 4; iBase++) {
                if (assignment.presentBase[iBase] && !presentBase[iBase]) {
                    addThis = false;
                }
            }
            if (addThis) {
                localAssignmentVector.add(assignment);
            }
        }
        smellTest(false);
        return siteActive;
    }
    
    int getBestGuessNPresent(){
        return bestGuessNPresent;
    }
    
    int[] getTotReads() {
        return totReads;
    }

        
    boolean smellTest(boolean printExceptions) {
        if (true) {
            return false;
        }
        double estProb = 0.0;
        double actProb = 0.0;
        boolean smellBad = false;
        for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
            if (totReads[iTimePoint] * totStrand[iTimePoint][0]*totStrand[iTimePoint][1] > 0) {
                for (int iBase = 0; iBase < 4; iBase++) {
                    if (reads[iTimePoint][iBase] > 0) {
                        estProb += this.gammaCalc.logGamma(reads[iTimePoint][iBase] + 0.5) 
                                - this.gammaCalc.logGamma(reads[iTimePoint][iBase] + 1.0) - 0.5723649;
                        actProb += -reads[iTimePoint][iBase] * 0.6931472 + this.gammaCalc.logGamma(reads[iTimePoint][iBase]+1)
                                - this.gammaCalc.logGamma(1.0 + 0.5 * totReads[iTimePoint] * strandReads[iTimePoint][0][iBase] / totStrand[iTimePoint][0]) 
                                - this.gammaCalc.logGamma(1.0 + 0.5 * totReads[iTimePoint] * strandReads[iTimePoint][1][iBase] / totStrand[iTimePoint][1]);
                    }
                }
            }
        }
        smellBad = (actProb - estProb < -11.512); // 11.512
        if (printExceptions && smellBad) {
            System.out.println("Rejected site\t" + iSite + "\t" + (actProb - estProb) + "\t" + smellBad);
        }
        return smellBad;
    }    

    
}
