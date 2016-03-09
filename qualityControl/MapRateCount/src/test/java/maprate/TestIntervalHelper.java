package maprate;

import htsjdk.samtools.cram.index.CramIndex;
import junit.framework.TestCase;

import java.util.List;
import java.util.Map;

/**
 * Created by Administrator on 2016/3/9.
 */
public class TestIntervalHelper extends TestCase {
    public void testBuild(){
        String filepath = "C:\\Users\\Administrator\\Desktop\\CCP.20131001.designed.bed";
        IntervalHelper bedHelper = new IntervalHelper(filepath);
        Map<String, List<IntervalHelper.Interval>> map = bedHelper.getIntervalMap();
        for(Map.Entry<String, List<IntervalHelper.Interval>> entry: map.entrySet()){
            System.out.println(entry.getKey());
            for(IntervalHelper.Interval interval : entry.getValue()){
                //System.out.println(interval);
            }
        }
    }

    public void testGetOverlap(){
        String filepath = "C:\\Users\\Administrator\\Desktop\\CCP.20131001.designed.bed";
        IntervalHelper bedHelper = new IntervalHelper(filepath);
        assertEquals(bedHelper.getOverlapLength("chr1", 10000, 10200), 0);
        assertEquals(bedHelper.getOverlapLength("chr1", 2488064, 2488070), 3);
        assertEquals(bedHelper.getOverlapLength("chr1", 2488068, 2488068), 1);
        assertEquals(bedHelper.getOverlapLength("chr1", 2489273, 2489772), 2);
        assertEquals(bedHelper.getOverlapLength("chr1", 2489273, 2489792), 22);
    }

}
