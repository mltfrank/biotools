package maprate;

import htsjdk.samtools.SAMRecord;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Created by Administrator on 2016/3/9.
 */
public class ReadCounter {
    // read statistic
    int unmappedReadCount = 0;
    int totalReadCount = 0;
    int secondaryOrSupplementaryReadCount = 0;

    // overlap statistic
    int totalOverlapLength = 0;
    int totalReadLength = 0;
    int overlapReadCount = 0;

    /**
     * Count read statistic
     * @param record
     */
    public void countRead(SAMRecord record){
        totalReadCount ++;
        if(record.isSecondaryOrSupplementary())
            secondaryOrSupplementaryReadCount ++;
        if(record.getReadUnmappedFlag())
            unmappedReadCount ++;
    }

    /**
     * Count overlap statistic
     * @param overlapLength read overlapped length
     * @param readLength read length
     */
    public void countOverlap(int overlapLength, int readLength){
        if(overlapLength > 0){
            totalOverlapLength += overlapLength;
            overlapReadCount ++;
        }
        totalReadLength += readLength;
    }

    private String getPercentageString(float val){
        val *= 100;
        return val + "%";
    }

    public void generateReport(String reportFilePath, boolean countOverlap){
        StringBuffer sb = new StringBuffer();
        sb.append("Total reads: ");
        sb.append(totalReadCount);
        sb.append("\nMapped read count: ");
        sb.append(totalReadCount-unmappedReadCount-secondaryOrSupplementaryReadCount);
        sb.append("  ");
        sb.append(getPercentageString((float)(totalReadCount-unmappedReadCount-secondaryOrSupplementaryReadCount)
                / (float)totalReadCount));
        sb.append("\nUnmapped read count: ");
        sb.append(unmappedReadCount);
        sb.append("  ");
        sb.append(getPercentageString((float)unmappedReadCount / (float)totalReadCount));
        sb.append("\nSecondary or supplementary read: ");
        sb.append(secondaryOrSupplementaryReadCount);
        sb.append("  ");
        sb.append(getPercentageString((float)secondaryOrSupplementaryReadCount / (float)totalReadCount));

        // if need overlap report
        if(countOverlap){
            sb.append("\nTotal length of reads: ");
            sb.append(totalReadLength);
            sb.append("\nRead has overlaps: ");
            sb.append(overlapReadCount);
            sb.append("  ");
            sb.append(getPercentageString((float)overlapReadCount / (float)totalReadCount));
            sb.append("\nTotal overlap length in reads: ");
            sb.append(totalOverlapLength);
            sb.append("  ");
            sb.append(getPercentageString((float)totalOverlapLength / (float)totalReadLength));
        }
        try {
            FileWriter writer = new FileWriter(new File(reportFilePath));
            writer.write(sb.toString());
            writer.flush();
            writer.close();
        }catch (IOException e){
            e.printStackTrace();
        }
    }
}
