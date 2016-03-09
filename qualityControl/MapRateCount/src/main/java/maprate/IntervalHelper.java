package maprate;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;

import java.io.*;
import java.util.*;

/**
 * Read interval file or bed file.
 * Author: Wang Bingchen
 * Date: 2016/3/9.
 */
public class IntervalHelper {
    private final Log log = Log.getInstance(IntervalHelper.class);

    /**
     * A hash map for quick search interval.
     * This map use a contig name defined in reference's dict file ('1', '2', '3' in b37 reference,
     * or 'chr1', 'chr2' in hg19 reference) as key and map key to a list of intervals on this contig.
     * When given a range(start coordinate and end coordinate of a read) to find if the range has overlap
     * with one or more interval, this map can quickly get the list of intervals on this contig, and scan
     * list linearly.
     */
    Map<String, List<Interval>> intervalMap = null;

    /**
     * Build interval helper
     * @param intervalFilePath interval file path
     */
    public IntervalHelper(String intervalFilePath){
        if (!validFileType(intervalFilePath)) {
            throw new RuntimeException("Don't support file format: \'"+intervalFilePath+"\'");
        }
        this.intervalMap = readIntervalFile(intervalFilePath);
        if(this.intervalMap == null){
            throw new RuntimeException("Fail to read intervals");
        }
    }

    /**
     * Valid interval file type.
     * Only support .interval or .bed file in current version.
     * @param intervalFilePath interval file path
     * @return if file format can be read in this version.
     */
    private static boolean validFileType(String intervalFilePath){
        return (intervalFilePath.endsWith(".bed") || intervalFilePath.endsWith(".interval"));
    }

    /**
     * Read interval file to parse the interval, store it into a map.
     * The interval map is described above, as variable 'intervalMap'.
     * Read file as text(just as .bed or .interval file). Parse each line into a interval.
     * The rule to parse is as follow:
     *      The first column is the contig name.
     *      The second column is start coordinate.
     *      The third column is end coordinate.
     * If start coordinate is larger than end coordinate, swap them.
     * If a line is not fit the rule, just jump over the line.(In case of describe line and state line).
     * @param intervalFilePath
     * @return
     */
    private Map<String, List<Interval>> readIntervalFile(String intervalFilePath){
        Map<String, List<Interval>> newIntervalMap = new HashMap<String, List<Interval>>();
        try {
            BufferedReader fileReader = new BufferedReader(new FileReader(new File(intervalFilePath)));
            String line = fileReader.readLine();
            // traverse each line in file and parse to interval
            while(line != null){
                if(line.startsWith("@"))
                    continue;
                try{
                    String[] lineSplit = line.split("\t");
                    String contigName = lineSplit[0];
                    int startCord = Integer.parseInt(lineSplit[1]);
                    int endCord = Integer.parseInt(lineSplit[2]);

                    // put into map, group by chromosome name
                    List<Interval> chromList = newIntervalMap.get(contigName);
                    if(chromList == null){
                        chromList = new ArrayList<Interval>();
                        newIntervalMap.put(contigName, chromList);
                    }
                    Interval newInterval = (startCord < endCord) ?
                            new Interval(startCord, endCord) : new Interval(endCord, startCord);
                    chromList.add(newInterval);
                }catch (Exception e){
                    log.warn("Jump over line \'"+line+"\'");
                }
                line = fileReader.readLine();
            }
        } catch (FileNotFoundException e){
            log.error("File \'"+intervalFilePath+"\' is not found.");
            return null;
        } catch (IOException e){
            log.error("Fail to read file \'"+intervalFilePath+"\'");
            return null;
        }

        // sort interval by start
        for(List<Interval> list: newIntervalMap.values()){
            Collections.sort(list, new IntervalComparator());
        }

        return newIntervalMap;
    }

    public Map<String, List<Interval>> getIntervalMap(){
        return this.intervalMap;
    }

    public int getOverlapLength(SAMRecord record){
        String refName = record.getReferenceName();
        int readStart = record.getAlignmentStart();
        int readEnd = record.getAlignmentEnd();
        return (readStart < readEnd) ? getOverlapLength(refName, readStart, readEnd)
                : getOverlapLength(refName, readEnd, readStart);
    }

    /**
     * Judge if a read have overlap with some panel intervals.
     * A read may have overlap with more than one panels on same chromosome,
     * but can't have overlap with panels on different chromosome.
     *
     * Note: We assume the panel don't have each other. If have, the result length may be more than expected.
     *
     * @param chromName contig name same with the definition interval file.
     * @param start read alignment start
     * @param end read alignment end
     * @return overlap length in total.
     */
    public int getOverlapLength(String chromName, int start, int end){
        List<Interval> list = this.intervalMap.get(chromName);
        int overlapLength = 0;
        for(Interval interval : list){
            if(interval.start > end){ // As interval is sorted by start, if a panel's start is larger than end of a read,
                // the following interval will surely have no overlap with given read.
                return overlapLength;
            }
            overlapLength += interval.getOverlap(start, end);
        }
        return overlapLength;
    }

    /**
     * Comparator for interval sort.
     * Sorted by interval's start coordinate.
     */
    class IntervalComparator implements Comparator<Interval>{
        public int compare(Interval inter1, Interval inter2){
            return inter1.start - inter2.start;
        }
    }

    public class Interval{
        private int start;
        private int end;
        public Interval(int start, int end){
            this.start = start;
            this.end = end;
        }

        public int getOverlap(int newStart, int newEnd){
            int maxStart = (newStart>start) ? newStart : start;
            int minEnd = (newEnd<end) ? newEnd:end;
            if(maxStart <= minEnd)
                return minEnd - maxStart + 1;
            else
                return 0;
        }

        @Override
        public String toString(){
            return "start:" + start + "   end:"+end;
        }
    }
}
