/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package cluster_rg;


import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.apache.commons.math3.optim.*;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import picocli.CommandLine;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import pal.math.ConjugateDirectionSearch;
import pal.math.OrthogonalHints;
import pal.math.ConjugateGradientSearch;
import flanagan.math.Minimisation;
import flanagan.math.MinimisationFunction;




public class Cluster_RG {
    
    static PrintWriter hapWriter = null;                    // PrintWriter for haplotypes
    static PrintWriter logWriter = null;                    // PrintWriter for log files
    static PrintWriter logLikeWriter = null;                 // PrintWriter for likelihood terms
    
    double best = 1.0E20;
    double[] bestPoint = null;
    
    double[] lb_alpha = {0.95, 0.00001};
    double[] ub_alpha = {0.9999, 0.2};
    
    public static void main(String[] args) {
        Cluster_RG m = new Cluster_RG();
        m.run(args);
    }

    private void run(String[] args) {
        Options options = new Options();
        CommandLine cmd = new CommandLine(options);

        try {
            cmd.parse(args);
            if (cmd.isUsageHelpRequested()) {
                cmd.usage(System.err);
            } else if (cmd.isVersionHelpRequested()) {
                cmd.printVersionHelp(System.err);
            } else {
                long startTime = System.currentTimeMillis();
                validateOptions(options);
                
                initialiseFiles(options, args);
                
                GammaCalc gammaCalc = GammaCalc.get(options.gammaCache);

                // fraction of sites to use when optimising alpha parameters
                Constants.USE_FRAC[0] = options.alpha_frac;
                Constants.USE_FRAC[1] = options.alpha_frac;

                long fileSeed = options.randomSeed;

                List<Cluster> clusters = new ArrayList<>();
                for (int i = 0; i < options.countFile.length; i++) {
                    File initialFreqFile;
                    if (options.initialFreqFile == null) {
                        initialFreqFile = null;
                    } else {
                        initialFreqFile = options.initialFreqFile[i];
                    }
                    Cluster cluster = new Cluster(options.countFile[i],
                            initialFreqFile,
                            options.haplotypes[i],
                            options.initialAlphaParams,
                            options.optimiser,
                            gammaCalc,
                            fileSeed++,
                            options.verbose);
                    clusters.add(cluster);
                }

                // Optimise
                optimise(clusters, options);

                long endTime = System.currentTimeMillis();
                System.out.printf("Main: Execution time = %.2fs\n", (endTime - startTime) / 1000.0);
                logWriter.printf("Main: Execution time = %.2fs\n", (endTime - startTime) / 1000.0);
                
                try{
                    logWriter.close();
                    if (options.printHaplotypes) {
                        hapWriter.close();
                    }
                    if (options.printLikelihoods) {
                        logLikeWriter.close();
                    }                    
                } catch(Exception ex) {
                    System.exit(1);
                }
                
            }
            

        } catch (CommandLine.ParameterException ex) {
            System.err.println(ex.getMessage());
            logWriter.println(ex.getMessage());
            ex.getCommandLine().usage(System.err);
        } catch (Exception ex) {
            throw new CommandLine.ExecutionException(cmd, "Error", ex);
        }

    }

    private void validateOptions(Options options) {
        if (options.countFile.length != options.haplotypes.length) {
            String msg = String.format("You have %d files but %d haplotype numbers.\n", options.countFile.length, options.haplotypes.length);
            throw new RuntimeException(msg);
        }

        if (options.initialAlphaParams == null) {
            options.initialAlphaParams = new double[]{Constants.DEFAULT_ALPHA_0, Constants.DEFAULT_ALPHA_1};
        }
    }

    private void optimise(List<Cluster> clusters, Options options) {
        final ExecutorService threadPool = Executors.newFixedThreadPool(options.threads);
        ConvergenceChecker<PointValuePair> convergenceChecker = new SimpleValueChecker(-1, options.tol);
        final double[] currentAlphaParams = options.initialAlphaParams;
 
        if (options.process) {            
            System.out.println(Arrays.toString(currentAlphaParams));
            for (Cluster cluster : clusters) {
                double[][] hapParams = cluster.getCurrentHapParams();
                for (int iDim = 0; iDim < hapParams.length; iDim++) {
                    System.out.println(Arrays.toString(hapParams[iDim]));
                }
            }
        } else {
            System.out.format("%s%.8f %.8f\n",
                    "Starting pre-optimisation of hap frequencies, a = ", 
                    currentAlphaParams[0], currentAlphaParams[1]);
            logWriter.format("%s%.8f %.8f\n",
                    "Starting pre-optimisation of hap frequencies, a = ", 
                    currentAlphaParams[0], currentAlphaParams[1]);
            List<Future<Double>> futures = new ArrayList<>();
            for (final Cluster cluster : clusters) {
                Future<Double> future = threadPool.submit(() -> cluster.run());
                futures.add(future);
            }
            List<Double> output = Cluster_RG.getFutureResults(futures);
            double totalLogLikelihood = output.stream().mapToDouble(Double::doubleValue).sum();
            if (options.fixAlpha) {
                System.out.println("Finished haplo-freq optimsation, Total loglikelihood = " + totalLogLikelihood);
                logWriter.println("Finished haplo-freq optimsation, Total loglikelihood = " + totalLogLikelihood);
            } else {
                System.out.println("Finished pre-optimsation, Total loglikelihood = " + totalLogLikelihood);
                logWriter.println("Finished pre-optimsation, Total loglikelihood = " + totalLogLikelihood);
                double[] tempAlpha = optimiseAlpha(clusters, currentAlphaParams, threadPool);
                currentAlphaParams[0] = tempAlpha[0];
                currentAlphaParams[1] = tempAlpha[1];
            }
        }
 
        threadPool.shutdown();

        System.out.println("\nMain: Converged.");
        System.out.println("\n\n========================= RESULTS =========================");
        logWriter.println("\nMain: Converged.");
        logWriter.println("\n\n========================= RESULTS =========================");

        double finalLnl = 0;
        for (Cluster cluster : clusters) {
            System.out.println();
            logWriter.println();
            double lnlContrib = cluster.printResults(options.printHaplotypes, options.printLikelihoods);
            finalLnl += lnlContrib;
        }

        System.out.printf("\nMain: Final total likelihood = %.7f\n", finalLnl);
        logWriter.printf("\nMain: Final total likelihood = %.7f\n", finalLnl);
    }

    private double[] optimiseAlpha(List<Cluster> clusters, double[] startAlpha, ExecutorService threadPool) {
        if (true) {
            MinimisationFunction clusterAlphaOptimise = new OptimiseAlphaFunction(clusters, threadPool);
            double[] point = Arrays.copyOf(startAlpha, 2);
            Minimisation min = new Minimisation();
            double[] step = {0.001, 0.001};
            min.nelderMead(clusterAlphaOptimise, point, step, 0.1);
            min.addConstraint(0, -1, lb_alpha[0]);
            min.addConstraint(0, 1, ub_alpha[0]);
            min.addConstraint(1, -1, lb_alpha[1]);
            min.addConstraint(1, 1, ub_alpha[1]);
            min.getMinimum();
            return point;
        } else if (false) {
            MultivariateOptimizer optimize = new BOBYQAOptimizer(2 * 2, 0.0001, 1.0E-7);
            org.apache.commons.math3.analysis.MultivariateFunction optimiseAlphaFunction 
                    = new OptimiseAlphaFunction(clusters, threadPool);
            OptimizationData[] parm = new OptimizationData[]{       // Set up optimisation data
                    new InitialGuess(Arrays.copyOf(startAlpha,2)),
                    new MaxEval(1000000),
                    GoalType.MINIMIZE,
                    new ObjectiveFunction(optimiseAlphaFunction),
                    new SimpleBounds(lb_alpha, ub_alpha)};
            double[] optPoint = optimize.optimize(parm).getPoint();  // Optimise
            return optPoint;
        } else {
            pal.math.MultivariateFunction clusterAlphaOptimise = new OptimiseAlphaFunction(clusters, threadPool);      
            ConjugateDirectionSearch cds = new ConjugateDirectionSearch();
            cds.step = 0.001;
            double[] point = Arrays.copyOf(startAlpha, 2);
            cds.optimize(clusterAlphaOptimise, point, 0.01, 0.1E-6);
            return point;
        }

    }

    private class OptimiseAlphaFunction implements  org.apache.commons.math3.analysis.MultivariateFunction, pal.math.MultivariateFunction, MinimisationFunction {
        final List<Cluster> clusters;
        final ExecutorService threadPool;
        private OptimiseAlphaFunction(final List<Cluster> clusters, ExecutorService threadPool) {
            this.clusters = clusters;
            this.threadPool = threadPool;
        }

        @Override
        public double evaluate(double[] point) { 
                return value(point);
        }
        @Override
        public double getLowerBound(int n) { 
                return lb_alpha[n];
        }    
        @Override
        public double getUpperBound(int n) { 
                return ub_alpha[n];
        } 
        
        @Override
        public int getNumArguments() { 
                return 2;
        }    
         
        @Override
        public OrthogonalHints getOrthogonalHints() { 
                return null;
        }    
        
        public double function(double[] point) {
            for (int iPoint = 0; iPoint < point.length; iPoint++) {
                if (point[iPoint] < getLowerBound(iPoint) || point[iPoint] > getUpperBound(iPoint)) {
                    return 1.0E20;
                }
            }
            return value(point);
        }
        
        public double value(double[] point) {
            List<Future<Double>> futures = new ArrayList<>();

            for (final Cluster cluster : clusters) {
                // optimise each cluster haplotypes independently (no synchronisation req)
                Future<Double> future = threadPool.submit(() -> cluster.optimiseAlpha(1, point));
                futures.add(future);
            }
            List<Double> output = Cluster_RG.getFutureResults(futures);
            double tot = output.stream().mapToDouble(Double::doubleValue).sum();
            if (tot < best) {
                best = tot;
                bestPoint = Arrays.copyOf(point, 2);
            }
            System.out.println("Optimising alpha = " + Arrays.toString(point) + "\tTotal loglikelihood = " 
                    + tot + "\tBest loglikelihood = " + best);

            return tot;
        }
    }
    
    void initialiseFiles(Options options, String[] args) {                
        try{
            String logFileName = options.prefix + ".log";
            File logFile = new File(logFileName);
            if (logFile.exists()){
                logFile.delete();
            }
            logWriter = new PrintWriter(logFile);    

            if (options.printHaplotypes) {
                String hapFileName = options.prefix + "Haplo.fasta";
                File hapFile = new File(hapFileName);
                if (hapFile.exists()){
                    hapFile.delete();
                }
                hapWriter = new PrintWriter(hapFile);                         
            }

            if (options.printLikelihoods) {
                String logLikeFileName = options.prefix + ".lld";
                File logLikeFile = new File(logLikeFileName);
                if (logLikeFile.exists()){
                    logLikeFile.delete();
                }
                logLikeWriter = new PrintWriter(logLikeFile);                         
            }                

        } catch (IOException ioe) {
            System.out.println("Error -- " + ioe.toString());
            logWriter.println("Error -- " + ioe.toString());
            System.exit(1);
        }

        System.out.printf("Main: arguments = %s\n", String.join(" ", args));
        logWriter.printf("Main: arguments = %s\n", String.join(" ", args));
        // Setup
        System.out.printf("Main: seed = %d\n", options.randomSeed);
        logWriter.printf("Main: seed = %d\n", options.randomSeed);
        System.out.println("Optimisation method :" + options.optimiser);
        logWriter.println("Optimisation method :" + options.optimiser);
    }
    
    
    private static <T> List<T> getFutureResults(List<Future<T>> futures) {
        List<T> results = new ArrayList<>();

        for (Future<T> f : futures) {
            try {
                results.add(f.get());
            } catch (Exception e) {
                e.printStackTrace();
                throw new RuntimeException(e);
            }
        }
        return results;
    }

}