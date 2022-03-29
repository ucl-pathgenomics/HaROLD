package cluster_rg;

import flanagan.math.Minimisation;
import flanagan.math.MinimisationFunction;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import org.apache.commons.math3.exception.MathIllegalStateException;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import pal.math.ConjugateDirectionSearch;
import pal.math.OrthogonalHints;

/**
 *
 * @author rgoldst
 */
public class Cluster {

    private final int nHaplo; // Number of haplotypes, revised based on command line argument
    private int nTimePoints = 0;  // Number of timepoints, revised based on data
    private int[] nActive;
    
    private ArrayList<Assignment> assignmentVector = new ArrayList<>();  // Vector of all possible assignments
    private int[] nAssignDiffBases = null; // Number of assignments with a given number of bases
    private DataSet dataSet;  // Class for holding and manipulating sequence data

    private Random random;
    private boolean verbose; // Print lots of intermediate results

    private int maxIter = Constants.MAX_ITER; // Maximum rounds of optimisation
    private double minImprovement = Constants.DEFAULT_TOL;  // Minimum improvement necessary to continue optimisation
    private double[] initialAlphaParams;
    
    private String optimiser;
    
    private double finalLogLikelihood = 0.0;
    private double currentLogLikelihood;

    private final String name;
    private int nEvaluations = 0;
    
    private boolean firstPassage = true;
    
    private double[][] currentHapParams;
    private double[] currentAlphaParams;
    private HashMap<double[], double[][]> previousHapFreqHash = new HashMap<>();
    private double[] scale ={1.0, 1.0};
    private double[] previousAlphaParams = new double[2];
    
    /**
    * Reads in data and initialises
    */  
    Cluster(File countFilesFile, File initialFreqFile, int nHaplo, double[] initialAlpha,
            String optimiser, GammaCalc gammaCalc, long randomSeed, boolean verbose) {
        nAssignDiffBases = new int[5];
        this.name = countFilesFile.getName();
        this.optimiser = optimiser;
        System.out.println(this.name + ": " + countFilesFile.getAbsolutePath());
        Cluster_RG.logWriter.println(this.name + ": " + countFilesFile.getAbsolutePath());
        this.random = new Random(randomSeed);
        this.verbose = verbose;

        this.initialAlphaParams = initialAlpha;

        this.nHaplo = nHaplo;  // Update number of haplotypes
        System.out.printf("%s: haplotypes = %d\n", this.name, this.nHaplo);
        Cluster_RG.logWriter.printf("%s: haplotypes = %d\n", this.name, this.nHaplo);

        constructAssignments(gammaCalc);  // Construct possible assignments of bases to haplotypes
        dataSet = new DataSet(countFilesFile, initialFreqFile, nHaplo, assignmentVector, nAssignDiffBases,
                gammaCalc, random, verbose); // Construct dataset
        nTimePoints = dataSet.getNTimePoints();  // Number of time points in dataset
        nActive = dataSet.getNActive();
        this.currentHapParams = dataSet.getInitialHapParams();  // Start with initial nearly equal haplotype frequencies
        this.currentAlphaParams = Arrays.copyOf(initialAlphaParams, 2);   // Initial values for alpha parameters alpha0 and alphaE
        System.out.println(Arrays.toString(nActive));
        System.out.printf("%s: timepoints = %d\n", this.name, this.nTimePoints);
        System.out.printf("%s: sites = %d\n", this.name, dataSet.getSiteCount());
        Cluster_RG.logWriter.println(Arrays.toString(nActive));
        Cluster_RG.logWriter.printf("%s: timepoints = %d\n", this.name, this.nTimePoints);
        Cluster_RG.logWriter.printf("%s: sites = %d\n", this.name, dataSet.getSiteCount());
        previousHapFreqHash.put(this.initialAlphaParams, currentHapParams); 
    }

    
    double[][] getCurrentHapParams() {
        return currentHapParams;
    }
    

    /**
    * Find best assignments and haplotype frequencies
    */   
    double run() {
        int iIter = 0;
        dataSet.setParams(currentHapParams, currentAlphaParams);
        double step1_current_lnl = dataSet.value();
        double step1_previous_lnl = Double.NEGATIVE_INFINITY;

        while (iIter < maxIter) {

            if (Math.abs(step1_current_lnl - step1_previous_lnl) < minImprovement) {
                break;
            }
            nEvaluations++;
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) { 

                if (nActive[iTimePoint] == 0) {
                    currentHapParams[iTimePoint][0] = 1.0;
                    dataSet.setOptType(1, iTimePoint, iIter); 
                    Cluster_RG.logWriter.println("Fixed\t" + name + "\t" + iTimePoint + "\t" + iIter + "\t" + Arrays.toString(dataSet.computePiHap(iTimePoint, currentHapParams[iTimePoint])));
                } else if (nActive[iTimePoint] == 1) {   // Simple single parameter optimisation for each time point
                    double optSinglePoint;      //  Hapltype frequency parameter
                    dataSet.setOptType(1, iTimePoint, iIter);    // Tell dataSet what timePoint is being optimised
                    Cluster_RG.logWriter.println("Initial\t" + name + "\t" + iTimePoint + "\t" + iIter + "\t" + Arrays.toString(dataSet.computePiHap(iTimePoint, currentHapParams[iTimePoint])));
                    optSinglePoint = fmin(0., 1.0, 1.0E-6);    // Find best value within range and tolerance
                    Cluster_RG.logWriter.println("Optimised\t" + name + "\t" + iTimePoint + "\t" + iIter + "\t" + Arrays.toString(dataSet.computePiHap(iTimePoint, currentHapParams[iTimePoint])));
                    if (verbose) {
                        Cluster_RG.logWriter.println("Optimum piParams\t" + name + "\t" + iTimePoint + "\t" + iIter + "\t" + optSinglePoint);  // Output optimum
                    }
                    currentHapParams[iTimePoint][0] = optSinglePoint; // Update current HapParameters
                } else if (nActive[iTimePoint] > 1) {   // Multidimensional parameter optimisation for each time point
                    dataSet.setOptType(1, iTimePoint, iIter);   // Tell dataSet what timePoint is being optimised
                    if (optimiser.equals("BOBYQ")) {
                        if (firstPassage) {
                            System.out.println("In Cluster: BOBYQ");
                            firstPassage = false;
                        }
                        MultivariateOptimizer optimize = new BOBYQAOptimizer(2 * nActive[iTimePoint], 0.01, 1.0E-6);
                        org.apache.commons.math3.analysis.MultivariateFunction optimiseHaploFrequency 
                                = new OptimiseHaploFrequency();
                        double[] lb_haplo = new double[nActive[iTimePoint]];
                        double[] ub_haplo = new double[nActive[iTimePoint]];
                        Arrays.fill(ub_haplo, 1.0);
                        Cluster_RG.logWriter.println("Initial\t" + name + "\t" + iTimePoint + "\t" + iIter + "\t" + 
                                Arrays.toString(dataSet.computePiHap(iTimePoint, currentHapParams[iTimePoint])));
                        OptimizationData[] parm = new OptimizationData[]{       // Set up optimisation data
                                new InitialGuess(Arrays.copyOf(currentHapParams[iTimePoint],nActive[iTimePoint])),
                                new MaxEval(1000000),
                                GoalType.MINIMIZE,
                                new ObjectiveFunction(optimiseHaploFrequency),
                                new SimpleBounds(lb_haplo, ub_haplo)};
                        double[] optPoint = null;
                        try{
                            optPoint = optimize.optimize(parm).getPoint();  // Optimise
                        } catch (MathIllegalStateException mise) {
                            optimiser = "conjugate";
                            break;
                        }
                        if (true || verbose) {
                            Cluster_RG.logWriter.println("Optimum piParams\t" + name + "\t" + iTimePoint + "\t" + 
                                    iIter + "\t" + Arrays.toString(optPoint));
                        }
                        for (int iHaplo = 0; iHaplo < nActive[iTimePoint]; iHaplo++) {
                            currentHapParams[iTimePoint][iHaplo] = optPoint[iHaplo]; // Update current parameters
                        }
                    } else if (optimiser.equals("simplex")) {
                         if (firstPassage) {
                            System.out.println("In Cluster: Simplex");
                            firstPassage = false;
                        }
                        Minimisation min = new Minimisation();
                        MinimisationFunction optimiseHaploFrequency = new OptimiseHaploFrequency();
                        double[] point = Arrays.copyOf(currentHapParams[iTimePoint],nActive[iTimePoint]);
                        double[] step = new double[nActive[iTimePoint]];
                        Arrays.fill(step, 0.1);
  
                        min.nelderMead(optimiseHaploFrequency, point, 0.1, 10000);
                        for (int iDim = 0; iDim < nActive[iTimePoint]; iDim++) {
                            min.addConstraint(iDim, -1, 0.0);
                            min.addConstraint(iDim, 1, 1.0);
                        }
                        min.getMinimum();
                        for (int iHaplo = 0; iHaplo < nActive[iTimePoint]; iHaplo++) {
                            currentHapParams[iTimePoint][iHaplo] = point[iHaplo]; // Update current parameters
                        }
                    } else if (optimiser.equals("conjugate")) {  
                        if (firstPassage) {
                            System.out.println("In Cluster: conjugate");
                            firstPassage = false;
                        }
                        ConjugateDirectionSearch cds = new ConjugateDirectionSearch();
                        double[] point = Arrays.copyOf(currentHapParams[iTimePoint],nActive[iTimePoint]);
                        pal.math.MultivariateFunction optimiseHaploFrequency 
                                = new OptimiseHaploFrequency(iTimePoint);
                        cds.optimize(optimiseHaploFrequency, point, 0.1, 0.1);
                        for (int iHaplo = 0; iHaplo < nActive[iTimePoint]; iHaplo++) {
                            currentHapParams[iTimePoint][iHaplo] = point[iHaplo]; // Update current parameters
                        }
                    }
                    Cluster_RG.logWriter.println("Optimised\t" + name + "\t" + iTimePoint + "\t" + iIter 
                            + "\t" + Arrays.toString(dataSet.computePiHap(iTimePoint, currentHapParams[iTimePoint]))
                            + "\t" + nEvaluations + "\t" + dataSet.nEvaluations);
                }
                
                
            }
            
            dataSet.setOptType(2, 0, 0);
            updatePriors();

            step1_previous_lnl = step1_current_lnl;
            step1_current_lnl = dataSet.value();  // Find best set of assignments
            System.out.println(name + "\tIteration " + iIter + "\tCurrent loglikelihood " + step1_current_lnl
                    + "\tImprovement " + (step1_current_lnl - step1_previous_lnl));
            Cluster_RG.logWriter.println(name + "\tIteration " + iIter + "\tCurrent loglikelihood " + step1_current_lnl
                    + "\tImprovement " + (step1_current_lnl - step1_previous_lnl));
            for (int iTimePoint = 0; iTimePoint < nTimePoints; iTimePoint++) {
                System.out.print(iTimePoint);
                Cluster_RG.logWriter.print(iTimePoint);
                double[] hapValues = dataSet.computePiHap(iTimePoint, currentHapParams[iTimePoint]);
                for (int iHaplo = 0; iHaplo < nHaplo; iHaplo++) {
                    System.out.format("\t%.5f", hapValues[iHaplo]);
                    Cluster_RG.logWriter.format("\t%.5f", hapValues[iHaplo]);
                }
                System.out.println();
                Cluster_RG.logWriter.println();
            }
            iIter++;

        }
        if (previousHapFreqHash.size() > 1000) {
            previousHapFreqHash.clear();
        }
        previousHapFreqHash.put(currentAlphaParams, currentHapParams);
        Cluster_RG.logWriter.printf("%s: haplotype frequencies lnl = %.5f\n", this.name, step1_current_lnl);
        return step1_current_lnl;
    }

    private class OptimiseHaploFrequency implements org.apache.commons.math3.analysis.MultivariateFunction, pal.math.MultivariateFunction, MinimisationFunction {
        
        int iTimePoint = 0;
        
        OptimiseHaploFrequency(int iTimePoint) {
            this.iTimePoint = iTimePoint;
        }
        
        OptimiseHaploFrequency(){
        }
        
        public double evaluate(double[] point) { 
                return dataSet.value(point);
        }

        public double getLowerBound(int n) { 
                return 0;
        }    

        public double getUpperBound(int n) { 
                return 1;
        } 
        
        public int getNumArguments() { 
                return nActive[iTimePoint];
        }    

        public OrthogonalHints getOrthogonalHints() { 
                return null;
        }    
        
        public double function(double[] point) {
            return dataSet.value(point);
        }
        
        public double value(double[] point) {
            return dataSet.value(point);
        }
    }


    
    void updatePriors() {
        this.dataSet.updatePriors();
    }

    double printResults(boolean printHaplotypes, boolean printLikelihoods) {
        return printResults(currentHapParams, currentAlphaParams, printHaplotypes, printLikelihoods);
    }
    
    
    double printResults(double[][] hapParams, double[] alphaParams, boolean printHaplotypes, boolean printLikelihoods) {
        System.out.println("Alpha level evaluations: " + nEvaluations + "\tHaplo level evaluations: " + dataSet.nEvaluations);
        System.out.println(Arrays.toString(alphaParams));
        Cluster_RG.logWriter.println(Arrays.toString(alphaParams));
        for (int i = 0; i < hapParams.length; i++) {
            System.out.println(i + "\t" + Arrays.toString(hapParams[i]));
            Cluster_RG.logWriter.println(i + "\t" + Arrays.toString(hapParams[i]));
        }
        dataSet.setOptType(2, 0, 0);
        dataSet.setParams(hapParams, alphaParams);
        finalLogLikelihood = dataSet.value();  // Find best set of assignments and calculate loglikelihood
        System.out.printf("-------------------- %s --------------------\n", this.name);
        Cluster_RG.logWriter.printf("-------------------- %s --------------------\n", this.name);
        dataSet.printResults(printHaplotypes, printLikelihoods);
        return finalLogLikelihood;
    }

    double optimiseAlpha(int iIter, double[] alphaParams) {
        this.currentAlphaParams = alphaParams;
        this.dataSet.setOptType(0, 0, iIter);   // Instruct dataSet to optimise alpha0 and alphaE
        if (previousHapFreqHash.size() > 1) {
            scale[0] = 0.9 * Math.abs(alphaParams[0] - previousAlphaParams[0]);
            scale[1] = 0.9 * Math.abs(alphaParams[1] - previousAlphaParams[1]);
            double bestDistance = 1.0E8;
            double distance = 0.0;
            for (double[] oldAlphas : previousHapFreqHash.keySet()) {
                distance = (alphaParams[0]-oldAlphas[0]) * (alphaParams[0]-oldAlphas[0]) / (scale[0]*scale[0]) 
                        + (alphaParams[1]-oldAlphas[1]) * (alphaParams[1]-oldAlphas[1]) / (scale[1]*scale[1]);
                if (distance < bestDistance) {
                    bestDistance = distance;
                    this.currentHapParams = previousHapFreqHash.get(oldAlphas);
                }
            }
        }
        previousAlphaParams = Arrays.copyOf(alphaParams, 2);
        double val = this.run();
        return -val;
    }

    
    /**
    * Constructs vector of all possible assignments
    */       
    private void constructAssignments(GammaCalc gammaCalc) {
        int nAssignments = pow(Constants.MAX_BASES, nHaplo);  // Theoretical exhaustive number of possible assignments
        for (int iAssign = 0; iAssign < nAssignments; iAssign++) {  // Loop over all possible assignments
            Assignment newAssignment = new Assignment(iAssign, nHaplo, gammaCalc, verbose);
            assignmentVector.add(newAssignment);
            nAssignDiffBases[newAssignment.nPresent]++;
        }
        if (false) {
            System.out.printf("%s: assignments = %d\n", name, assignmentVector.size());
            for (int i = 0; i < 5; i++) {
                System.out.println(i + "\t" + nAssignDiffBases[i]);
            }
        }
    }

    
    /**
    * Computes a^b
    */
    private int pow (int a, int b) {  // Computes powers
        if ( b == 0)     return 1;
        if ( b == 1)     return a;
        if (b%2 == 0)    return     pow ( a * a, b/2); //even a=(a^2)^b/2
        else             return a * pow ( a * a, (b-1)/2); //odd  a=a*(a^2)^b/2
    }
    
    /**
    * Brent algorithm for maximising a function in one dimension
    * Adapted from Apache Commons
    */    
    private double fmin (double a, double b, double tol) {
        double c,d,e,eps,xm,p,q,r,tol1,t2, u,v,w,fu,fv,fw,fx,x,tol3;

        c = .5*(3.0 - Math.sqrt(5.0));
        d = 0.0;
        eps = 1.2e-16;
        tol1 = eps + 1.0;
        eps = Math.sqrt(eps);

        v = a + c*(b-a);
        w = v;
        x = v;
        e = 0.0;
        fx=dataSet.value(x);
        fv = fx;
        fw = fx;
        tol3 = tol/3.0;

        xm = .5*(a + b);
        tol1 = eps*Math.abs(x) + tol3;
        t2 = 2.0*tol1;

        while (Math.abs(x-xm) > (t2 - .5*(b-a))) {
            p = q = r = 0.0;
            if (Math.abs(e) > tol1) {
                r = (x-w)*(fx-fv);
                q = (x-v)*(fx-fw);
                p = (x-v)*q - (x-w)*r;
                q = 2.0*(q-r);
                if (q > 0.0) {
                        p = -p;
                } else {
                        q = -q;
                }
                r = e;
                e = d;
            }

            if ((Math.abs(p) < Math.abs(.5*q*r)) && (p > q*(a-x)) && (p < q*(b-x))) {
                d = p/q;
                u = x+d;
                if (((u-a) < t2) || ((b-u) < t2)) {
                    d = tol1;
                    if (x >= xm) d = -d;
                }
            } else {
                if (x < xm) {
                    e = b-x;
                } else {
                    e = a-x;
                }
                d = c*e;
            }

            if (Math.abs(d) >= tol1) {
                u = x+d;
            } else {
                if (d > 0.0) {
                    u = x + tol1;
                } else {
                    u = x - tol1;
                }
            }
            fu = dataSet.value(u);

            if (fx <= fu) {
                if (u < x) {
                    a = u;
                } else {
                    b = u;
                }
            }
            if (fu <= fx) {
                if (u < x) {
                    b = x;
                } else {
                    a = x;
                }
                v = w;
                fv = fw;
                w = x;
                fw = fx;
                x = u;
                fx = fu;
                xm = .5*(a + b);
                tol1 = eps*Math.abs(x) + tol3;
                t2 = 2.0*tol1;
            } else {
                if ((fu <= fw) || (w == x)) {
                    v = w;
                    fv = fw;
                    w = u;
                    fw = fu;
                    xm = .5*(a + b);
                    tol1 = eps*Math.abs(x) + tol3;
                    t2 = 2.0*tol1;
                } else if ((fu > fv) && (v != x) && (v != w)) {
                    xm = .5*(a + b);
                    tol1 = eps*Math.abs(x) + tol3;
                    t2 = 2.0*tol1;
                } else {
                    v = u;
                    fv = fu;
                    xm = .5*(a + b);
                    tol1 = eps*Math.abs(x) + tol3;
                    t2 = 2.0*tol1;
                }
            }
        }
        return x;
    }

}
