/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package refineHaplotypes;

import java.io.File;
import java.util.Arrays;

/**
 *
 * @author rgoldst
 */
public class ReadArguments {
    private String baseFreqFileName = "";
    private String hapAlignmentFileName = "";
    private String hapFreqFileName = "";
    private String readsFileName = "";
    private String refSeqFileName = "";
    private String tag = "";
    private File baseFreqFile;
    private File hapAlignmentFile;
    private File hapFreqFile;
    private File readsFile;
    private File refSeqFile;
    private String[] arguments;
    private double minReadDepth = -1.;
    private double minReads = 20.;
    private double errorRate = 0.002;                                     // Expected error rate
    private boolean iterate = false;
    private int maxIter = 10;
    private int maxHaplo = 10;
    private int maxRecombine = 20;
    private boolean expand = false;
       
    ReadArguments(String[] args) {
        this.arguments = Arrays.copyOf(args, args.length);
        int iArg = 0;
        while (iArg < args.length) {
            if (args[iArg].startsWith("-")){
                if (args[iArg].equals("-bam") || args[iArg].equals("-sam" )) {
                    this.readsFileName = args[iArg+1]; 
                    iArg+=2;
                } else if (args[iArg].equals("-baseFreq")) {
                    this.baseFreqFileName = args[iArg+1];
                    iArg+=2; 
                } else if (args[iArg].equals("-hapFreq")) {
                    this.hapFreqFileName = args[iArg+1]; 
                    iArg+=2;
                } else if (args[iArg].equals("-hapSeq")) {
                    this.hapAlignmentFileName = args[iArg+1]; 
                    iArg+=2;
                } else if (args[iArg].equals("-refSeq")) {
                    this.refSeqFileName = args[iArg+1]; 
                    iArg+=2;
                } else if (args[iArg].equals("-tag")) {
                    this.tag = args[iArg+1];
                    iArg+=2;
                } else if (args[iArg].equals("-minReadDepth")) {
                    this.minReadDepth = Double.parseDouble(args[iArg+1]); 
                    iArg+=2;
                } else if (args[iArg].equals("-minReads")) {
                    this.minReads = Double.parseDouble(args[iArg+1]); 
                    iArg+=2;
                } else if (args[iArg].equals("-errorRate")) {
                    this.errorRate = Double.parseDouble(args[iArg+1]); 
                    iArg+=2;
                } else if (args[iArg].equals("-iterate")) {
                    this.iterate = true; 
                    iArg++;
                } else if (args[iArg].equals("-maxIterate")) {
                    this.maxIter = Integer.parseInt(args[iArg+1]); 
                    iArg+=2;
                } else if (args[iArg].equals("-maxHaplo")) {
                    this.maxHaplo = Integer.parseInt(args[iArg+1]); 
                    iArg+=2;
                } else if (args[iArg].equals("-Expand")) {
                    this.expand = true; 
                    iArg++;
                } else if (args[iArg].equals("-maxRecombine")) {
                    this.maxRecombine = Integer.parseInt(args[iArg+1]); 
                    iArg+=2;
                }
            }
        }
        if (maxIter < 1) {
            iterate = false;
        }
        boolean notOK = false;
        if (readsFileName.length() == 0) {
            System.out.println("Must supply BAM or SAM file name, \"-bam\" or \"-sam\"");
            notOK = true;
        }
        if (baseFreqFileName.length() == 0) {
            System.out.println("Must supply base frequency file name, \"-baseFreq\"");
            notOK = true;
        }
        if (hapAlignmentFileName.length() == 0) {
            System.out.println("Must supply hap alignment file name, \"-hapSeq\"");
            notOK = true;
        }
        if (refSeqFileName.length() == 0) {
            System.out.println("Must supply reference sequence file name, \"-refSeq\"");
            notOK = true;
        }
        if (tag.length() == 0) {
            System.out.println("Must supply output file tag, \"-tag\"");
            notOK = true;
        }
        if (!notOK) {
            baseFreqFile = new File(baseFreqFileName);
            hapAlignmentFile = new File(hapAlignmentFileName);
            hapFreqFile = new File(hapFreqFileName);
            readsFile = new File(readsFileName);
            refSeqFile = new File(refSeqFileName);     
            if (!baseFreqFile.exists()) {
                System.out.println(baseFreqFileName + " does not exist");
                notOK = true;
            }
            if (!hapAlignmentFile.exists()) {
                System.out.println(hapAlignmentFileName + " does not exist");
                notOK = true;
            }
            if (!readsFile.exists()) {
                System.out.println(readsFileName + " does not exist");
                notOK = true;
            }
            if (!refSeqFile.exists()) {
                System.out.println(refSeqFileName + " does not exist");
                notOK = true;
            }            
            
        }
        if (notOK) {
            System.exit(1);
        }
      
    }

    
    
    public File getBaseFreqFile() {
        return baseFreqFile;
    }

    
    public File getHapAlignmentFile() {
        return hapAlignmentFile;
    }

    
    public File getHapFreqFile() {
        return hapFreqFile;
    }

    
    public File getReadsFile() {
        return readsFile;
    }    
    
    
    public File getRefSeqFile() {
        return refSeqFile;
    }
    
    public String getTag() {
        return (tag);
    }

    
    public double getMinReadDepth() {
        return (minReadDepth);
    }
    
        
    public double getMinReads() {
        return (minReads);
    }

        
    public double getErrorRate() {
        return (errorRate);
    }
    
           
    public boolean getIterate() {
        return (iterate);
    }

        
    public int getMaxIter() {
        return (maxIter);
    }
        
        
    public int getMaxHaplo() {
        return (maxHaplo);
    }
    
        
    public int getMaxRecombine() {
        return (maxRecombine);
    }
    
        
    public boolean getExpand() {
        return (expand);
    }
    

    
    public void printStuff() {
        System.out.print("Arguments:");
        RefineHaplotypes.logWriter.print("Arguments:");
        for (int iArg = 0; iArg < arguments.length; iArg++) {
            System.out.print(" " + arguments[iArg]);
            RefineHaplotypes.logWriter.print(" " + arguments[iArg]);
        }
        System.out.println();
        RefineHaplotypes.logWriter.println();
        System.out.println("Minimum read depth: " + minReadDepth);
        System.out.println("Minimum reads: " + minReads);
        RefineHaplotypes.logWriter.println("Minimum read depth: " + minReadDepth);
        RefineHaplotypes.logWriter.println("Minimum reads: " + minReads);
        System.out.println("Error rate: " + errorRate);
        RefineHaplotypes.logWriter.println("Error rate: " + errorRate);
        System.out.println("Iterate: " + iterate);
        RefineHaplotypes.logWriter.println("Iterate: " + iterate);
        if (iterate) {
            System.out.println("Max haplotypes: " + maxHaplo);
            RefineHaplotypes.logWriter.println("Max haplotypes: " + maxHaplo);
            System.out.println("Max iterations: " + maxIter);
            RefineHaplotypes.logWriter.println("Max iterations: " + maxIter);
            System.out.println("Max recombine: " + maxRecombine);
            RefineHaplotypes.logWriter.println("Max recombine: " + maxRecombine);
            System.out.println("Expand: " + expand);
            RefineHaplotypes.logWriter.println("Expand: " + expand);
        }
        System.out.println();
        RefineHaplotypes.logWriter.println();
    }
    
    
    
}
