/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package cluster_rg;

/**
 *
 * @author rgoldst
 */
public class Fodder {
        
//    void setParams(int iStage, double alpha_e, double beta_e, double S, double F0, double[][] alpha){
//        alpha_e_loc = alpha_e;
//        beta_e_loc = beta_e;
//        S_loc = S;
//        F0_loc = F0;
//        alpha_loc = alpha;
//        if (iStage == 0) {
////            double[][] logLikeAssignPoly = new double[nTimePoints][nAssignments];
//            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
//                if (timePointHasData[iTimePoint]){
//                    for (int iAssignment = 1; iAssignment < nAssignments-1; iAssignment++) {
//                        double[] sumAlpha = new double[2];
//                        for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
//                            sumAlpha[assign[iAssignment][iHaplo]] += alpha_loc[iTimePoint][iHaplo];
//                        }
//                        logLikeAssignPoly[iTimePoint][iAssignment] = 
//                            logStrandChoose[iTimePoint][0] + logStrandChoose[iTimePoint][1]
//                            + Beta.logBeta(reads[iTimePoint][1]+sumAlpha[0], reads[iTimePoint][0]+sumAlpha[1])
//                            - Beta.logBeta(sumAlpha[0], sumAlpha[1]);
//                    }
//                }
//            }
//        }
//    }
//    
//    
//    double computeLogLikelihoodAlpha(double[] alpha, int iTimePoint) { // alpha[][] = iTP, iHaplo, iTP
//        if (currentAssignment == 0 || currentAssignment == nAssignments-1 || !timePointHasData[iTimePoint]) {
//            return 0.0;
//        }
//        double totalProb = 0.0;
//        double logProb = 0.0;           
//        double[] sumAlpha = new double[2];
//        for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
//            sumAlpha[assign[currentAssignment][iHaplo]] += alpha[iHaplo];
//        }
//        logProb = Math.log(probStrand[iTimePoint][0]) + Math.log(probStrand[iTimePoint][1])
//            + logStrandChoose[iTimePoint][0] + logStrandChoose[iTimePoint][1]
//            + Beta.logBeta(reads[iTimePoint][1]+sumAlpha[0], reads[iTimePoint][0]+sumAlpha[1])
//            - Beta.logBeta(sumAlpha[0], sumAlpha[1]);
//        return logProb;
//    }
//
//        
//    double computeLogLikelihood(double alpha_e, double beta_e, double S, double F0, boolean printHaplo) { // alpha[][] = iHaplo, iTP
//        double[] logLikeAssign = new double[nAssignments];  // log likelihood of the site for each assignment
//        double totalProb = 0.0;
//        assignmentProbs = new double[nAssignments];
//        double[][] likeAssignSNoS = new double[nTimePoints][2];
//        for (int iAssignment = 0; iAssignment < nAssignments; iAssignment++) {
//            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
//                if (!timePointHasData[iTimePoint]) {
//                    likeAssignSNoS[iTimePoint][0] = 1.0;
//                    likeAssignSNoS[iTimePoint][1] = 0.0;
//                }
//                if (timePointHasData[iTimePoint]) {
//                    double[][][] probStrandSNoS = new double[nAssignments][2][2];
//                    
//                    if (iAssignment == 0) {
//                        for (int iStrand = 0; iStrand < 2; iStrand++) {
//                            double logProbTerm = logStrandChoose[iTimePoint][iStrand]
//                            + Beta.logBeta(strandReads[iTimePoint][iStrand][1]+alpha_e, strandReads[iTimePoint][iStrand][0]+beta_e)
//                            - Beta.logBeta(alpha_e, beta_e);
//                            probStrandSNoS[iAssignment][iStrand][0] = (1.0-S)*Math.exp(logProbTerm);
//                            probStrandSNoS[iAssignment][iStrand][1] = S/(strandReads[iTimePoint][iStrand][0]+strandReads[iTimePoint][iStrand][1]+1.0);
//                            probStrand[iTimePoint][iStrand] = probStrandSNoS[iAssignment][iStrand][0] 
//                                    + probStrandSNoS[iAssignment][iStrand][1];
//                        }
//                        logLikeAssign[iAssignment] += Math.log(probStrand[iTimePoint][0]) + Math.log(probStrand[iTimePoint][1]);
//                        likeAssignSNoS[iTimePoint][0] += probStrandSNoS[iAssignment][0][0]*probStrandSNoS[iAssignment][1][0];
//                        likeAssignSNoS[iTimePoint][1] += probStrand[iTimePoint][0]*probStrand[iTimePoint][1]
//                                - probStrandSNoS[iAssignment][0][0]*probStrandSNoS[iAssignment][1][0];
//                    } else if (iAssignment == nAssignments - 1) {
//                        for (int iStrand = 0; iStrand < 2; iStrand++) {
//                            double logProbTerm = logStrandChoose[iTimePoint][iStrand]
//                                + Beta.logBeta(strandReads[iTimePoint][iStrand][0]+alpha_e, strandReads[iTimePoint][iStrand][1]+beta_e)
//                                - Beta.logBeta(alpha_e, beta_e);
//                            probStrandSNoS[iAssignment][iStrand][0] = (1.0 - S) * Math.exp(logProbTerm);
//                            probStrandSNoS[iAssignment][iStrand][1] = S/(strandReads[iTimePoint][iStrand][0] + strandReads[iTimePoint][iStrand][1]+1.0);
//                            probStrand[iTimePoint][iStrand] = probStrandSNoS[iAssignment][iStrand][0] 
//                                    + probStrandSNoS[iAssignment][iStrand][1];
//                        }
//                        logLikeAssign[iAssignment] += Math.log(probStrand[iTimePoint][0]) + Math.log(probStrand[iTimePoint][1]);
//                        likeAssignSNoS[iTimePoint][0] += probStrandSNoS[iAssignment][0][0]*probStrandSNoS[iAssignment][1][0];
//                        likeAssignSNoS[iTimePoint][1] += probStrand[iTimePoint][0]*probStrand[iTimePoint][1]
//                                - probStrandSNoS[iAssignment][0][0]*probStrandSNoS[iAssignment][1][0];
//
//                    } else {
//                        for (int iStrand = 0; iStrand < 2; iStrand++) {
//                            double logProbTerm = 
//                                Beta.logBeta(alpha_e, strandReads[iTimePoint][iStrand][0]+strandReads[iTimePoint][iStrand][1] + beta_e) 
//                                    - Beta.logBeta(alpha_e, beta_e);
//                             probStrand[iTimePoint][iStrand] = (1.0-S)*Math.exp(logProbTerm)
//                                     + S/(strandReads[iTimePoint][iStrand][0]+strandReads[iTimePoint][iStrand][1]+1.0);
//                        }
//                        logLikeAssign[iAssignment] += Math.log(probStrand[iTimePoint][0]) + Math.log(probStrand[iTimePoint][1])
//                                + logLikeAssignPoly[iTimePoint][iAssignment];
//                        likeAssignSNoS[iTimePoint][0] += probStrand[iTimePoint][0]*probStrand[iTimePoint][1]*Math.exp(logLikeAssignPoly[iTimePoint][iAssignment]);
//                    }   // end of elses
//                }  // end of conditional on timepoint existing
//            }   // end of loop over timepoints
// 
//            if (iAssignment == 0 || iAssignment == nAssignments - 1) {
//                assignmentProbs[iAssignment] = (1.0 - F0) * 0.5 * Math.exp(logLikeAssign[iAssignment]);
//            } else {
//                assignmentProbs[iAssignment] += (F0 / (nAssignments - 2.0)) * Math.exp(logLikeAssign[iAssignment]);
//            }
//            totalProb += assignmentProbs[iAssignment];
//        } // end of loop over assignments
//        currentAssignment = iSelect(assignmentProbs);
//        if (printHaplo) {
//            double maxNoSProb = 1.0;
//            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
//                maxNoSProb *= likeAssignSNoS[iTimePoint][0]/(likeAssignSNoS[iTimePoint][0]+likeAssignSNoS[iTimePoint][1]);
//            }
//            System.out.print(iSite);
//            for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
//                double[] sumProb = new double[2];
//                for (int iAssignment = 0; iAssignment < nAssignments; iAssignment++) {
//                    sumProb[assign[iAssignment][iHaplo]] += assignmentProbs[iAssignment];
//                }
//            }
//            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
//                if (iTimePoint == 0) {
//                    System.out.print("\t");
//                } else {
//                    System.out.print(" ");
//                }
//                System.out.print(Arrays.toString(reads[iTimePoint]));
//            }
//            System.out.println();
//        }
//        return Math.log(totalProb);
//    }
//    
//    int iSelect(double[] probs) {   // Selects integer from a specified array of probabilities
//        int nPoss = probs.length;
//        double summ = 0.0;
//        for (int iPoss = 0; iPoss < nPoss; iPoss++) {
//            summ += probs[iPoss];
//        }
//        for (int iPoss = 0; iPoss < nPoss; iPoss++) {
//            probs[iPoss] /= summ;
//        }
//        double randomNumber = rag.Cluster.random.nextDouble();
//        for (int iPoss = 0; iPoss < nPoss-1; iPoss++) {
//            if (probs[iPoss] > randomNumber) {
//                return iPoss;
//            } else {
//                randomNumber -= probs[iPoss];
//            }
//        }
//        return (nPoss-1);
//    }

}
