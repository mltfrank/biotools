package maprate;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;

/**
 * Main class to run map rate counter
 * Generate report file.
 * Author: Wang Bingchen
 * Date: 2016/3/8.
 */
public class MapRateCounter {
    private final Log log = Log.getInstance(MapRateCounter.class);

    /** @argument input sam file */
    private String inputFile = null;
    /** @argument report file */
    private String reportFile = null;
    /** @argument interval file to statistic total base in panel. */
    private String intervalFile = null;

    private void parseArguments(String[] args){
        for(String arg : args){
            if(arg.startsWith("INPUT=")){
                this.inputFile = arg.substring(6);
            }
            else if(arg.startsWith("OUTPUT=")){
                this.reportFile = arg.substring(7);
            }
            else if(arg.startsWith("INTERVAL=")){
                this.intervalFile = arg.substring(9);
                log.info("Found interval file \'"+this.intervalFile+"\', program will count total overlap in report");
            }
            else{
                throw new RuntimeException("Unsupported arguments : \'"+arg+"\'");
            }
        }
        // check arguments
        if(inputFile == null){
            printUsage();
            throw new RuntimeException("Missing input file in arguments");
        }
        if(reportFile == null){
            printUsage();
            throw new RuntimeException("Missing report file in arguments");
        }
    }

    private void printUsage(){
        System.out.println("Usage:");
        System.out.println("java -jar MapRateCounter INPUT=<inputFile> OUTPUT=<outputFile> [INTERVAL=<intervalFile>]");
        System.out.println("  INPUT=<inputFile> required, sam or bam file path to be count");
        System.out.println("  OUTPUT=<outputFile> required, report file path");
        System.out.println("  INTERVAL=<intervalFile> optional, interval file path");
    }

    private void run(String[] args){
        if(args.length == 0) {
            printUsage();
            return;
        }
        // parse and check arguments
        parseArguments(args);

        // counter
        IntervalHelper intervalHelper = null;
        ReadCounter counter = new ReadCounter();

        // if need to count overlap, read interval file.
        boolean countOverlap = (this.intervalFile != null);
        if(countOverlap){
            log.info("Reading intervals from file: "+intervalFile);
            intervalHelper = new IntervalHelper(this.intervalFile);
            log.info("Reading intervals from file: "+intervalFile);
        }

        log.info("Reading records from file: "+inputFile);
        IOHelper.SamHeaderAndIterator headerAndIterator = IOHelper.openInput(inputFile);
        final SAMFileHeader header = headerAndIterator.header;
        final CloseableIterator<SAMRecord> iterator = headerAndIterator.iterator;

        log.info("Start traversing reads from file..");
        while(iterator.hasNext()){
            SAMRecord record = iterator.next();
            counter.countRead(record);
            if(countOverlap)
                counter.countOverlap(intervalHelper.getOverlapLength(record), record.getReadLength());
        }
        iterator.close();
        log.info("Finish traversing reads");
        log.info("Generating report into file: "+reportFile);
        counter.generateReport(reportFile, countOverlap);
        log.info("Finish all task");
    }

    public static void main(String[] args){
        new MapRateCounter().run(args);
    }


}
