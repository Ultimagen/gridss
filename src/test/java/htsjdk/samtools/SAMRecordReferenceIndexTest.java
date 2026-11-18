package htsjdk.samtools;

import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;

public class SAMRecordReferenceIndexTest {
    private static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(Arrays.asList(
                new SAMSequenceRecord("chr1", 1000), new SAMSequenceRecord("chr2", 2000), new SAMSequenceRecord("chrUn_x", 30))));
        return header;
    }

    @Test
    public void indexSettersMatchNameSetters() {
        final SAMFileHeader header = header();
        for (int i = -1; i < header.getSequenceDictionary().size(); i++) {
            final String name = i < 0 ? SAMRecord.NO_ALIGNMENT_REFERENCE_NAME : header.getSequence(i).getSequenceName();
            final SAMRecord byIndex = new SAMRecord(header);
            byIndex.setReferenceIndex(i);
            byIndex.setMateReferenceIndex(i);
            final SAMRecord byName = new SAMRecord(header);
            byName.setReferenceName(new String(name));
            byName.setMateReferenceName(new String(name));
            Assert.assertEquals(byName.getReferenceIndex(), byIndex.getReferenceIndex());
            Assert.assertSame(byName.getReferenceName(), byIndex.getReferenceName());
            Assert.assertEquals(byName.getMateReferenceIndex(), byIndex.getMateReferenceIndex());
            Assert.assertSame(byName.getMateReferenceName(), byIndex.getMateReferenceName());
            Assert.assertEquals(byName.getSAMString(), byIndex.getSAMString());
        }
    }

    @Test
    public void indexSettersReplaceAPreviouslyReadIndex() {
        final SAMRecord r = new SAMRecord(header());
        r.setReferenceIndex(0);
        r.setMateReferenceIndex(0);
        Assert.assertEquals(0, r.getReferenceIndex().intValue());
        Assert.assertEquals(0, r.getMateReferenceIndex().intValue());
        r.setReferenceIndex(1);
        r.setMateReferenceIndex(2);
        Assert.assertEquals(1, r.getReferenceIndex().intValue());
        Assert.assertEquals(2, r.getMateReferenceIndex().intValue());
    }

    @Test
    public void unknownIndexThrows() {
        final SAMRecord r = new SAMRecord(header());
        Assert.assertThrows(IllegalArgumentException.class, () -> r.setReferenceIndex(3));
        Assert.assertThrows(IllegalArgumentException.class, () -> r.setMateReferenceIndex(3));
    }

    @Test
    public void noHeaderAcceptsOnlyNoAlignment() {
        final SAMRecord r = new SAMRecord(null);
        r.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
        r.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
        Assert.assertEquals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME, r.getReferenceName());
        Assert.assertEquals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME, r.getMateReferenceName());
        Assert.assertThrows(IllegalStateException.class, () -> r.setReferenceIndex(0));
        Assert.assertThrows(IllegalStateException.class, () -> r.setMateReferenceIndex(0));
    }
}
