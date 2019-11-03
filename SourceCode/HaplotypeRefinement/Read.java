package refineHaplotypes;

import htsjdk.samtools.SAMRecord;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 * Class for containing information about a single read
 * @author rgoldst
 */
public class Read{
    private int iStart;                         // start position
    private int iLength;                        // length
    private int[] sequence;                     // sequence of bases: [lSite] (e.g. local site)
    private byte[] baseQual;
    private boolean[] siteExists;
    private int[] limits = new int[2];          // start and end position
    private boolean negativeStrand;             // direction of strand
    private HashMap<Integer, Integer> sigSiteHash = new HashMap<>();
    private int nSigSites = 0;
    private int nCopies = 1;
    private String sigSiteTag = "";

    /**
     * Constructor for Read
     * Sets parameters describing read including sequence
     * Default sequence value is 4 for no information e.g. '-'
     * 
     * @param samRecord samRecord to be stored
     */
    public Read(SAMRecord samRecord) {
        this.iStart = samRecord.getAlignmentStart();
        this.iLength = samRecord.getAlignmentEnd() - samRecord.getAlignmentStart();
        this.negativeStrand = samRecord.getReadNegativeStrandFlag();
        limits[0] = iStart;
        limits[1] = iStart+iLength;
        this.sequence = new int[iLength];
        this.baseQual = new byte[iLength];
        this.siteExists = new boolean[iLength];
        Arrays.fill(this.sequence, 4);
        Arrays.fill(this.baseQual, (byte) 0);
        Arrays.fill(this.siteExists, false);
        for (int iSite = iStart; iSite < iStart + iLength; iSite++) {
            int readPosition = samRecord.getReadPositionAtReferencePosition(iSite)-1;
            if (readPosition >= 0 && 
                    (samRecord.getBaseQualities().length == 0 || 
                            samRecord.getBaseQualities()[readPosition] > RefineHaplotypes.minBaseQual)) {
                this.siteExists[iSite-iStart] = true;             
                this.sequence[iSite-iStart] = RefineHaplotypes.getDNA(samRecord.getReadString().charAt(readPosition));
            }
        }
    }

    /**
     * Returns hashmap of bases in read at significant sites
     * 
     * @return HashMap of bases at significant sites
     */
    public HashMap<Integer, Integer> getSigSiteHash() {
        return sigSiteHash;
    }
    
    /**
     * Set the significant sites for this read and return a String containing list of these sites
     * 
     * @param sigSiteList List of 'significant' variable sites
     * @param basePresent List of acceptable bases for that position
     * @return String representing significant sites for this read
     */
    public void setSignificantSites(ArrayList<Integer> sigSiteList) {
        for (int iSite : sigSiteList) {
            int lSite = iSite - iStart;
            // if a site is significant, add to sigSiteHash and tag
            if ((iSite >= iStart) && (iSite < (iStart+iLength))
                    && siteExists[lSite] && (sequence[lSite] < 4)) {
                    sigSiteHash.put(iSite, sequence[lSite]);
                    if (sigSiteTag.length() > 0) {
                        sigSiteTag = sigSiteTag + "," + iSite + ":" + sequence[lSite];
                    } else {
                        sigSiteTag = sigSiteTag + "" + iSite + ":" + sequence[lSite];
                    }
            }
        }
        nSigSites = sigSiteHash.size();
    }
    
    public String getSigSiteTag() {
        return sigSiteTag;
    }
    
    public int getNCopies() {
        return nCopies;
    }
    
    public void addCopies() {
        nCopies++;
    }

    
    /**
     * Get the sequence of this read
     * @return sequence in integer format
     */
    public int[] getSequence() {
        return sequence;
    }

        
    /**
     * Get the sequence of this read
     * @return sequence in integer format
     */
    public boolean[] getSiteExists() {
        return siteExists;
    }
    
    /**
     * Get limits on read
     * @return start and stop locations
     */
    public int[] getLimits() {
        return limits;
    }
    
    /**
     * 
     * @param iSite site to be queried (local index)
     * @return sequence at this site
     */
    public int getBase(int iSite) {
        return sequence[iSite-iStart];
    }
    
    /**
     * Query strand direction
     * @return whether read is negative strand
     */
    public boolean getNegativeStrand() {
        return negativeStrand;
    }

}    
