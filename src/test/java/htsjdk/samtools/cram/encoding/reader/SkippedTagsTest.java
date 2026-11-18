package htsjdk.samtools.cram.encoding.reader;

import htsjdk.samtools.CRAMFileWriter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.cram.build.CramContainerIterator;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.cram.structure.CRAMEncodingStrategy;
import htsjdk.samtools.cram.structure.CompressionHeader;
import htsjdk.samtools.cram.structure.Container;
import htsjdk.samtools.cram.structure.DataSeries;
import htsjdk.samtools.cram.structure.EncodingDescriptor;
import htsjdk.samtools.cram.structure.EncodingID;
import htsjdk.samtools.cram.structure.ReadTag;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import org.junit.Assert;
import org.junit.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

public class SkippedTagsTest {
    private static final File REFERENCE = new File("src/test/resources/small.fa");

    private static CompressionHeader writeAndReadHeader() {
        final SAMFileHeader header = new SAMFileHeader(SAMSequenceDictionaryExtractor.extractDictionary(REFERENCE.toPath()));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        header.addReadGroup(new SAMReadGroupRecord("rg"));
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        final CRAMFileWriter writer = new CRAMFileWriter(new CRAMEncodingStrategy(), out, null, true, new ReferenceSource(REFERENCE), header, "test.cram");
        for (int i = 0; i < 20; i++) {
            final SAMRecord r = new SAMRecord(header);
            r.setReadName("r" + i);
            r.setReferenceIndex(0);
            r.setAlignmentStart(1 + 10 * i);
            r.setCigarString("10M");
            r.setReadString("ACGTACGTAC");
            r.setBaseQualityString("IIIIIIIIII");
            r.setAttribute("RG", "rg");
            r.setAttribute("t0", "IIIII" + i);
            r.setAttribute("tp", new byte[] { (byte) i, 0, 1 });
            r.setAttribute("NM", i);
            writer.addAlignment(r);
        }
        writer.close();
        final CramContainerIterator it = new CramContainerIterator(new ByteArrayInputStream(out.toByteArray()));
        while (it.hasNext()) {
            final Container container = it.next();
            if (container.getCompressionHeader() != null) {
                return container.getCompressionHeader();
            }
        }
        throw new IllegalStateException("no data container");
    }

    private static int tagId(final String name, final char type) {
        return ReadTag.name3BytesToInt(new byte[] { (byte) name.charAt(0), (byte) name.charAt(1), (byte) type });
    }

    @Test
    public void skips_only_requested_tags_with_own_blocks() {
        final CompressionHeader header = writeAndReadHeader();
        final SkippedTags skipped = new SkippedTags(header, new HashSet<>(Arrays.asList("t0", "tp", "zz")));
        Assert.assertEquals(new HashSet<>(Arrays.asList(tagId("t0", 'Z'), tagId("tp", 'B'))), skipped.tagIds);
        Assert.assertFalse(skipped.contentIds.isEmpty());
        for (final DataSeries dataSeries : DataSeries.values()) {
            final EncodingDescriptor d = header.getEncodingMap().getEncodingDescriptorForDataSeries(dataSeries);
            if (d != null && d.getEncodingID() == EncodingID.EXTERNAL) {
                Assert.assertFalse(dataSeries.name(), skipped.contentIds.contains(dataSeries.getExternalBlockContentId()));
            }
        }
    }

    @Test
    public void skips_nothing_by_default() {
        final SkippedTags skipped = new SkippedTags(writeAndReadHeader(), Collections.emptySet());
        Assert.assertTrue(skipped.tagIds.isEmpty());
        Assert.assertTrue(skipped.contentIds.isEmpty());
    }

    @Test
    public void keeps_a_tag_whose_block_another_tag_reads() {
        final CompressionHeader header = writeAndReadHeader();
        final EncodingDescriptor t0 = header.getTagEncodingMap().get(tagId("t0", 'Z'));
        header.getTagEncodingMap().put(tagId("xx", 'Z'), t0);
        final SkippedTags skipped = new SkippedTags(header, new HashSet<>(Arrays.asList("t0")));
        Assert.assertTrue(skipped.tagIds.isEmpty());
    }
}
