package makereadcount;

public class ReadArguments
{
    static String inputFileName;
    static String tag;
    static int maxSize;
    static byte minBaseQual;
    static int minMappingQual;
    static double minBalanceEValue;
    
    ReadArguments(final String[] args) {
        ReadArguments.inputFileName = args[args.length - 1];
        ReadArguments.tag = ReadArguments.inputFileName.replaceFirst(".*\\/", "").replaceFirst(".bam", "").replaceFirst(".sam", "");
        for (int iArg = 0; iArg < args.length - 1; ++iArg) {
            if (args[iArg].startsWith("-")) {
                if (args[iArg].equals("-N")) {
                    ReadArguments.maxSize = Integer.parseInt(args[iArg + 1]);
                }
                else if (args[iArg].equals("-minMappingQual")) {
                    ReadArguments.minMappingQual = Integer.parseInt(args[iArg + 1]);
                }
                else if (args[iArg].equals("-minBaseQual")) {
                    ReadArguments.minBaseQual = (byte)Integer.parseInt(args[iArg + 1]);
                }
                else if (args[iArg].equals("-minBalanceEValue")) {
                    ReadArguments.minBalanceEValue = Double.parseDouble(args[iArg + 1]);
                }
            }
        }
    }
    
    public static int getMinMappingQual() {
        return ReadArguments.minMappingQual;
    }
    
    public static double getMinBalanceEValue() {
        return ReadArguments.minBalanceEValue;
    }
    
    public static byte getMinBaseQual() {
        return ReadArguments.minBaseQual;
    }
    
    public static int getMaxSize() {
        return ReadArguments.maxSize;
    }
    
    public static String getInputFileName() {
        return ReadArguments.inputFileName;
    }
    
    public static String getLogFileName() {
        return invokedynamic(makeConcatWithConstants:(Ljava/lang/String;)Ljava/lang/String;, ReadArguments.tag);
    }
    
    public static String getOutputFileName() {
        return invokedynamic(makeConcatWithConstants:(Ljava/lang/String;)Ljava/lang/String;, ReadArguments.tag);
    }
    
    static {
        ReadArguments.maxSize = 300000;
        ReadArguments.minBaseQual = 30;
        ReadArguments.minMappingQual = 10;
        ReadArguments.minBalanceEValue = 0.1;
    }
}
