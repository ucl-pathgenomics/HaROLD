package cluster_rg;

import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

import java.io.File;

@Command(name = "richards-haplotype-model", footer = "Copyright (c) 2018 Richard A Goldstein", description = "", version = "1.0")
public class Options {
    
    @Option(names = {"-a", "--initial-alpha"}, arity = "2", description = "")
    double[] initialAlphaParams = new double[]{Constants.DEFAULT_ALPHA_0, Constants.DEFAULT_ALPHA_1}; 
    
    @Option(names = {"-A", "--fix-alpha"}, required = false, description = "Fix alpha parameters")
    boolean fixAlpha = false;    

    @Option(names = {"--alpha-frac"}, required = false, description = "Fraction of sites to use to optimise error parameters")
    double alpha_frac = 1.0;
    
    // TODO: -c and -n should be n-arity inputs - make arrays
    @Option(names = {"-c", "--count-file"}, arity = "1..*", required = true, description = "file containing list of count files")
    File[] countFile;
    
    // TODO: -c and -n should be n-arity inputs - make arrays
    @Option(names = {"-f", "--initial-freq-file"}, arity = "1..*", description = "optional file containing hap frequency values")
    File[] initialFreqFile = null;  

    @Option(names = {"-g", "--gamma-cache"}, description = "", hidden = true)
    int gammaCache = 1000;

    @Option(names = {"-h", "-?", "--help"}, usageHelp = true, description = "give this help list")
    protected boolean helpRequested;

    @Option(names = {"-H", "--printHaplotypes"}, required = false, description = "Print haplotypes")
    boolean printHaplotypes = false;    

    @Option(names = {"-L", "--printLikelihoods"}, required = false, description = "Print likelihoods")
    boolean printLikelihoods = false;        

    @Option(names = {"-n", "--haplotypes"}, arity = "1..*", required = true, description = "number of haplotypes")
    int[] haplotypes;

    @Option(names = {"-N", "--noOpt"}, required = false, description = "Process without optimising")
    boolean process = false;   

    @Option(names = {"-o", "--optimiser"}, required = false, description = "Optimiser for haplotype frequencies")
    String optimiser = "BOBYQ";      

    @Option(names = {"-p", "--prefix"}, required = false, description = "Results file prefix")
    String prefix = "Results";

    @Option(names = {"-s", "--seed"}, description = "")
    long randomSeed = System.currentTimeMillis();

    @Option(names = {"--threads"})
    int threads = 1;

    @Option(names = {"--tol"})
    double tol = Constants.DEFAULT_TOL;
               
    @Option(names = {"-v", "--verbose"}, description = "")
    boolean verbose = false;

    @Option(names = {"-V", "--version"}, versionHelp = true, description = "")
    boolean versionRequested;

}

