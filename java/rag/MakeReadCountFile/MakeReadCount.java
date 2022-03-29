package makereadcount;

import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import java.io.IOException;
import java.io.File;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.distribution.BinomialDistribution;
import java.util.Iterator;
import java.util.TreeSet;
import java.util.TreeMap;
import java.util.HashMap;
import java.util.Random;
import java.util.ArrayList;
import java.io.PrintWriter;

public class MakeReadCount
{
    static PrintWriter outputWriter;
    static PrintWriter logWriter;
    ArrayList<Read> readList;
    ReadArguments readArguments;
    char[] dna;
    static Random random;
    private static final HashMap<Character, Integer> DNA_CHAR_TO_INT_DNA_HASH;
    
    public static void main(final String[] args) {
        final MakeReadCount m = new MakeReadCount(args);
        m.run();
    }
    
    private MakeReadCount(final String[] args) {
        this.readList = new ArrayList<Read>();
        this.dna = new char[] { 'A', 'C', 'G', 'T', '-', 'X', 'a', 'c', 'g', 't', 'x' };
        if (args.length == 0) {
            System.out.println("Arguments:");
            System.out.println("\tInitial number of haplotypes");
            System.out.println("\tNumber of sites");
            System.out.println("\tAlignment containing reference sequence as first sequence");
            System.out.println("\tList of SAM or BAM files containing data");
            System.exit(0);
        }
        this.readArguments = new ReadArguments(args);
        this.openFiles();
        this.fillHash();
    }
    
    private void run() {
        this.readData();
        this.process();
        this.finish();
    }
    
    public void process() {
        final TreeMap<Integer, Integer> mapQualityHisto = new TreeMap<Integer, Integer>();
        final TreeMap<Integer, Integer> baseQualityHisto = new TreeMap<Integer, Integer>();
        final int[][] pUnevenStrands = new int[2][21];
        final int[][][] strandCount = new int[ReadArguments.getMaxSize()][2][5];
        final TreeSet<Integer> sitesThatExist = new TreeSet<Integer>();
        final int[] strand = new int[2];
        for (final Read read : this.readList) {
            int iStrand = 0;
            if (read.getNegativeStrand()) {
                iStrand = 1;
            }
            final int[] array = strand;
            final int n = iStrand;
            ++array[n];
            final int mapQuality = read.getMapQual();
            int mqp1 = 1;
            if (mapQualityHisto.containsKey(mapQuality)) {
                mqp1 = mapQualityHisto.get(mapQuality) + 1;
            }
            mapQualityHisto.put(mapQuality, mqp1);
            if (mapQuality >= ReadArguments.minMappingQual) {
                final byte[] baseQuality = read.getBaseQual();
                final int[] sequence = read.getSequence();
                final boolean[] siteExists = read.getSiteExists();
                final int[] limits = read.getLimits();
                for (int iSite = limits[0]; iSite < limits[1]; ++iSite) {
                    final int lSite = iSite - limits[0];
                    if (siteExists[lSite]) {
                        final int bq = baseQuality[lSite];
                        int bqp1 = 1;
                        if (baseQualityHisto.containsKey(bq)) {
                            bqp1 = baseQualityHisto.get(bq) + 1;
                        }
                        baseQualityHisto.put(bq, bqp1);
                        sitesThatExist.add(iSite);
                        final int[] array2 = strandCount[iSite][iStrand];
                        final int n2 = sequence[lSite];
                        ++array2[n2];
                    }
                }
            }
        }
        int totBasesObserved = 0;
        final TreeMap<Integer, Double> omittedPositions = new TreeMap<Integer, Double>();
        MakeReadCount.outputWriter.println("Position,r_depth,Consensus_res,A_f,A_r,C_f,C_r,G_f,G_r,T_f,T_r");
        for (int iSite2 = 1; iSite2 < ReadArguments.getMaxSize(); ++iSite2) {
            final int[] iCount = { strandCount[iSite2][0][0] + strandCount[iSite2][0][1] + strandCount[iSite2][0][2] + strandCount[iSite2][0][3], strandCount[iSite2][1][0] + strandCount[iSite2][1][1] + strandCount[iSite2][1][2] + strandCount[iSite2][1][3] };
            final int total = iCount[0] + iCount[1];
            totBasesObserved += total;
            if (total > 0) {
                final double[] pValues = this.checkStrandCount(strandCount[iSite2]);
                double e = 2.0 * sitesThatExist.size() * pValues[0];
                int log10e = -20;
                if (e > 1.0E-20) {
                    log10e = Math.round((float)Math.round(Math.log10(e)));
                }
                if (log10e <= 0) {
                    final int[] array3 = pUnevenStrands[0];
                    final int n3 = -log10e;
                    ++array3[n3];
                }
                e = sitesThatExist.size() * pValues[1];
                log10e = -20;
                if (e > 1.0E-20) {
                    log10e = Math.round((float)Math.round(Math.log10(e)));
                }
                if (log10e <= 0) {
                    final int[] array4 = pUnevenStrands[1];
                    final int n4 = -log10e;
                    ++array4[n4];
                }
                int iMax = 0;
                int maxCount = strandCount[iSite2][0][0] + strandCount[iSite2][1][0];
                if (strandCount[iSite2][0][1] + strandCount[iSite2][1][1] > maxCount) {
                    iMax = 1;
                    maxCount = strandCount[iSite2][0][1] + strandCount[iSite2][1][1];
                }
                if (strandCount[iSite2][0][2] + strandCount[iSite2][1][2] > maxCount) {
                    iMax = 2;
                    maxCount = strandCount[iSite2][0][2] + strandCount[iSite2][1][2];
                }
                if (strandCount[iSite2][0][3] + strandCount[iSite2][1][3] > maxCount) {
                    iMax = 3;
                    maxCount = strandCount[iSite2][0][3] + strandCount[iSite2][1][3];
                }
                int chooseStrand = 0;
                if (iCount[1] > iCount[0]) {
                    chooseStrand = 1;
                }
                if (e < ReadArguments.getMinBalanceEValue()) {
                    final double frac0 = (iCount[0] - strandCount[iSite2][0][iMax]) / (iCount[0] + 1.0E-10);
                    final double frac2 = (iCount[1] - strandCount[iSite2][1][iMax]) / (iCount[1] + 1.0E-10);
                    if (frac2 < frac0) {
                        chooseStrand = 1;
                    }
                }
                MakeReadCount.outputWriter.print(invokedynamic(makeConcatWithConstants:(IIC)Ljava/lang/String;, iSite2, iCount[chooseStrand], this.dna[iMax]));
                for (int iBase = 0; iBase < 4; ++iBase) {
                    MakeReadCount.outputWriter.print(invokedynamic(makeConcatWithConstants:(I)Ljava/lang/String;, strandCount[iSite2][chooseStrand][iBase]));
                }
                MakeReadCount.outputWriter.println();
            }
        }
        MakeReadCount.logWriter.println(invokedynamic(makeConcatWithConstants:(Ljava/lang/Object;)Ljava/lang/String;, sitesThatExist.last()));
        MakeReadCount.logWriter.println(invokedynamic(makeConcatWithConstants:(I)Ljava/lang/String;, sitesThatExist.size()));
        MakeReadCount.logWriter.format("Avg read depth: %.4f\n", 1.0 * totBasesObserved / (1.0 * sitesThatExist.last()));
        MakeReadCount.logWriter.format("Minimum mapping quality: %d\n", ReadArguments.getMinMappingQual());
        MakeReadCount.logWriter.format("Minimum base quality: %d\n", ReadArguments.getMinBaseQual());
        MakeReadCount.logWriter.format("Minimum balance E value: %.2g\n\n", ReadArguments.getMinBalanceEValue());
        MakeReadCount.logWriter.println("log10 E values for uneven strand distribution:");
        for (int iP = 0; iP < 20; ++iP) {
            MakeReadCount.logWriter.println(invokedynamic(makeConcatWithConstants:(III)Ljava/lang/String;, -iP, pUnevenStrands[0][iP], pUnevenStrands[1][iP]));
        }
        MakeReadCount.logWriter.println(invokedynamic(makeConcatWithConstants:(II)Ljava/lang/String;, pUnevenStrands[0][20], pUnevenStrands[1][20]));
        MakeReadCount.logWriter.println();
        MakeReadCount.logWriter.println("Mapping quality histogram:");
        for (final int mq : mapQualityHisto.keySet()) {
            MakeReadCount.logWriter.println(invokedynamic(makeConcatWithConstants:(ILjava/lang/Object;)Ljava/lang/String;, mq, mapQualityHisto.get(mq)));
        }
        MakeReadCount.logWriter.println();
        MakeReadCount.logWriter.println("Base quality histogram:");
        for (final int bq2 : baseQualityHisto.keySet()) {
            MakeReadCount.logWriter.println(invokedynamic(makeConcatWithConstants:(ILjava/lang/Object;)Ljava/lang/String;, bq2, baseQualityHisto.get(bq2)));
        }
        MakeReadCount.logWriter.println();
        MakeReadCount.logWriter.println(invokedynamic(makeConcatWithConstants:(II)Ljava/lang/String;, strand[0], strand[1]));
        MakeReadCount.logWriter.println();
        if (omittedPositions.size() > 0) {
            for (final int iSite3 : omittedPositions.keySet()) {
                MakeReadCount.logWriter.format("%d\t%.2g\n", iSite3, omittedPositions.get(iSite3));
            }
        }
    }
    
    double[] checkStrandCount(final int[][] sCount) {
        final double[] pValues = new double[2];
        final int[] iCount = { sCount[0][0] + sCount[0][1] + sCount[0][2] + sCount[0][3], sCount[1][0] + sCount[1][1] + sCount[1][2] + sCount[1][3] };
        final int total = iCount[0] + iCount[1];
        int iMin = 0;
        if (iCount[1] < iCount[0]) {
            iMin = 1;
        }
        final BinomialDistribution bd = new BinomialDistribution(total, 0.5);
        pValues[0] = bd.cumulativeProbability(iCount[iMin]);
        pValues[1] = 100.0;
        if (iCount[0] > 2 && iCount[1] > 2) {
            for (int iBase = 0; iBase < 4; ++iBase) {
                final int sumBase = sCount[0][iBase] + sCount[1][iBase];
                if (sumBase > 0 && sumBase < (total + 1) / 2) {
                    final double realTerm = CombinatoricsUtils.binomialCoefficientDouble(iCount[0], sCount[0][iBase]) * CombinatoricsUtils.binomialCoefficientDouble(iCount[1], sCount[1][iBase]);
                    final double[] p = new double[2];
                    for (int iTerm = 0; iTerm <= Math.min(sumBase, iCount[0]); ++iTerm) {
                        if (sumBase - iTerm <= iCount[1]) {
                            final double term = CombinatoricsUtils.binomialCoefficientDouble(iCount[0], iTerm) * CombinatoricsUtils.binomialCoefficientDouble(iCount[1], sumBase - iTerm);
                            if (term <= realTerm + 1.0E-5) {
                                final double[] array = p;
                                final int n = 0;
                                array[n] += term;
                            }
                            else {
                                final double[] array2 = p;
                                final int n2 = 1;
                                array2[n2] += term;
                            }
                        }
                    }
                    pValues[1] = Math.min(p[0] / (p[0] + p[1]), pValues[1]);
                }
            }
        }
        return pValues;
    }
    
    public void openFiles() {
        try {
            final File logFile = new File(ReadArguments.getLogFileName());
            if (logFile.exists()) {
                logFile.delete();
            }
            MakeReadCount.logWriter = new PrintWriter(logFile);
            final File outputFile = new File(ReadArguments.getOutputFileName());
            if (outputFile.exists()) {
                outputFile.delete();
            }
            MakeReadCount.outputWriter = new PrintWriter(outputFile);
        }
        catch (IOException ioe) {
            System.out.println(invokedynamic(makeConcatWithConstants:(Ljava/lang/String;)Ljava/lang/String;, ioe.toString()));
            MakeReadCount.logWriter.println(invokedynamic(makeConcatWithConstants:(Ljava/lang/String;)Ljava/lang/String;, ioe.toString()));
            System.exit(1);
        }
    }
    
    public void finish() {
        try {
            MakeReadCount.logWriter.close();
            MakeReadCount.outputWriter.close();
        }
        catch (Exception ex) {
            System.exit(1);
        }
    }
    
    public void readData() {
        try {
            final File inputFile = new File(ReadArguments.getInputFileName());
            this.readFromFile(inputFile);
        }
        catch (IOException e) {
            System.out.println("Error: File not found (IO error)");
            System.exit(1);
        }
    }
    
    private void readFromFile(final File inputSamOrBamFile) throws IOException {
        final SamReader reader = SamReaderFactory.makeDefault().open(inputSamOrBamFile);
        for (final SAMRecord samRecord : reader) {
            final Read newRead = new Read(samRecord);
            this.readList.add(newRead);
        }
    }
    
    private void fillHash() {
        MakeReadCount.DNA_CHAR_TO_INT_DNA_HASH.put('A', 0);
        MakeReadCount.DNA_CHAR_TO_INT_DNA_HASH.put('C', 1);
        MakeReadCount.DNA_CHAR_TO_INT_DNA_HASH.put('G', 2);
        MakeReadCount.DNA_CHAR_TO_INT_DNA_HASH.put('T', 3);
        MakeReadCount.DNA_CHAR_TO_INT_DNA_HASH.put('-', 4);
        MakeReadCount.DNA_CHAR_TO_INT_DNA_HASH.put('N', 4);
    }
    
    public static int getDNA(final char DNA) {
        return MakeReadCount.DNA_CHAR_TO_INT_DNA_HASH.get(DNA);
    }
    
    static {
        MakeReadCount.outputWriter = null;
        MakeReadCount.logWriter = null;
        MakeReadCount.random = new Random();
        DNA_CHAR_TO_INT_DNA_HASH = new HashMap<Character, Integer>();
    }
}
