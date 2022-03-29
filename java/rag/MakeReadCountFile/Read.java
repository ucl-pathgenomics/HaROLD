package makereadcount;

import java.util.Arrays;
import htsjdk.samtools.SAMRecord;

public class Read
{
    private int iStart;
    private int iLength;
    private int[] limits;
    private boolean negativeStrand;
    private byte[] baseQual;
    private int mapQual;
    private int[] sequence;
    private boolean[] siteExists;
    private boolean hasBaseQualities;
    SAMRecord samRecord;
    
    public Read(final SAMRecord samRecord) {
        this.limits = new int[2];
        this.samRecord = samRecord;
        this.iStart = samRecord.getAlignmentStart();
        this.iLength = samRecord.getAlignmentEnd() - samRecord.getAlignmentStart();
        this.negativeStrand = samRecord.getReadNegativeStrandFlag();
        this.mapQual = samRecord.getMappingQuality();
        this.limits[0] = this.iStart;
        this.limits[1] = this.iStart + this.iLength;
        this.sequence = new int[this.iLength];
        this.baseQual = new byte[this.iLength];
        this.siteExists = new boolean[this.iLength];
        Arrays.fill(this.sequence, 4);
        Arrays.fill(this.baseQual, (byte)0);
        Arrays.fill(this.siteExists, false);
        this.hasBaseQualities = false;
        if (samRecord.getBaseQualities().length > 0) {
            this.hasBaseQualities = true;
        }
        for (int iSite = this.iStart; iSite < this.iStart + this.iLength; ++iSite) {
            final int readPosition = samRecord.getReadPositionAtReferencePosition(iSite) - 1;
            if (readPosition >= 0) {
                this.siteExists[iSite - this.iStart] = true;
                if (this.hasBaseQualities) {
                    this.baseQual[iSite - this.iStart] = samRecord.getBaseQualities()[readPosition];
                }
                else {
                    this.baseQual[iSite - this.iStart] = ReadArguments.getMinBaseQual();
                }
                this.sequence[iSite - this.iStart] = MakeReadCount.getDNA(samRecord.getReadString().charAt(readPosition));
            }
        }
    }
    
    public int getSequence(final int iSite) {
        return this.sequence[iSite - this.iStart];
    }
    
    public int[] getSequence() {
        return this.sequence;
    }
    
    public int getMapQual() {
        return this.mapQual;
    }
    
    boolean[] getSiteExists() {
        return this.siteExists;
    }
    
    public byte[] getBaseQual() {
        return this.baseQual;
    }
    
    public boolean getNegativeStrand() {
        return this.negativeStrand;
    }
    
    public int[] getLimits() {
        return this.limits;
    }
    
    public int getStart() {
        return this.limits[0];
    }
}
