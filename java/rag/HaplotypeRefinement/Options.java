import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import java.io.File;



@Command(name = "richards-haplotype-model", footer = "Copyright (c) 2018 Richard A Goldstein", description = "", version = "1.0")
public class Options {
    
        
    
    double defaultCutOff = 6.5;
    private String defaultIdiosyncratic = "-1.0";
    private String defaultUnfoldedStates = "1.0E60";        
    private String defaultTemperature = "0.6";   
    

    @Option(names = {"--bam", "--BAM", "--sam", "--SAM"}, required = true, description = "Name of bam/sam file (req)")
    File readsFile = null;  
    
    @Option(names = {"--baseFreq"}, required = true, description = "Name of base frequency file (req)")
    File baseFreqFile = null;  
    
    @Option(names = {"--hapFreq"}, required = false, description = "Name of haplotype frequency file (null)")
    File hapFreqFile = new File("Null");  
     
    @Option(names = {"--hapAlignment", "--alignment"}, required = true, description = "Name of haplotype alignment file, fasta (req)")
    File hapAlignmentFile = null;  
     
    @Option(names = {"--refSequence", "--referenceSequence", "--reference"}, required = true, description = "Name of file containing reference sequence, fasta (req)")
    File refSeqFile = null;  
    
    @Option(names = {"-t", "--tag"}, required = true, description = "Tag for output files (req)")
    String tag;    

    @Option(names = {"-D", "--minReadDepth"}, required = false, description = "Minimum read depth (<0 indicates inactive) (-1.)")
    double minReadDepth = -1.0;
    
    @Option(names = {"-m", "--minReads"}, required = false, description = "Minimum number of reads (20.0)")
    double minReads = 20.0;

    @Option(names = {"--errorRate"}, required = false, description = "Error rate (0.002)")
    double errorRate = 0.002;
    
    @Option(names = {"-I", "--iterate"}, required = false, description = "Turn on iteration (false)")
    boolean iterate = false;
        
    @Option(names = {"--printReference"}, required = false, description = "Print reference sequence (false)")
    boolean printRef = false;
        
    @Option(names = {"--printIntermediate"}, required = false, description = "Print intermediate sequences (false)")
    boolean printIntermediate = false;
  
    @Option(names = {"--maxIterate"}, required = false, description = "Maximum number of iterations (10)")
    int maxIter = 10;
    
    @Option(names = {"--maxRecombine"}, required = false, description = "Maximum number of recombination attempts (20)")
    int maxRecombine = 20;
    
    @Option(names = {"--maxHaplo", "--maxHaplotypes"}, required = false, description = "Maximum number of haplotypes (10)")
    int maxHaplo = 10;
    
    @Option(names = {"--expand"}, required = false, description = "Consider splits that increase number of haplotypes (false)")
    boolean expand = false;
   
    @Option(names = {"--seed"}, required = false, description = "Random number seed (-1, not specified)")
    long randomSeed = -1;
    
    @Option(names = {"-h", "-?", "--help"}, usageHelp = true, description = "give this help list")
    protected boolean helpRequested;

    @Option(names = {"-V", "--version"}, versionHelp = true, description = "")
    boolean versionRequested;
    
    void checkOptions() {
        if (maxIter < 1) {
            iterate = false;
        }
        
        
    }

    
}
