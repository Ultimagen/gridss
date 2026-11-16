package htsjdk.samtools;

import htsjdk.samtools.cram.build.CramContainerIterator;
import htsjdk.samtools.cram.ref.CRAMReferenceSource;
import htsjdk.samtools.cram.structure.Container;
import htsjdk.samtools.cram.structure.Slice;
import htsjdk.samtools.cram.structure.block.BlockCompressionMethod;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * embedded_ref.cram: samtools 1.19.2 view -C, CRAM 3.0, embed_ref=1, seqs_per_slice=30, slices_per_container=2, of 330
 * synthetic reads on a synthetic order-1 Markov reference, so the embedded reference blocks are rANS order 1 as in
 * Ultima CRAMs. embedded_ref.bam: htslib's decode of it (decode_md=0).
 */
public class CRAMIteratorEmbeddedReferenceTest {
    private static final File CRAM = new File("src/test/resources/embedded_ref.cram");
    private static final File BAM = new File("src/test/resources/embedded_ref.bam");

    private static final CRAMReferenceSource NO_REFERENCE = new CRAMReferenceSource() {
        @Override
        public byte[] getReferenceBases(final SAMSequenceRecord sequenceRecord, final boolean tryNameVariants) {
            throw new AssertionError("embedded reference expected");
        }

        @Override
        public byte[] getReferenceBasesByRegion(final SAMSequenceRecord sequenceRecord, final int zeroBasedStart, final int requestedRegionLength) {
            throw new AssertionError("embedded reference expected");
        }
    };

    private static String fieldsAndSortedTags(final SAMRecord r) {
        final String[] f = r.getSAMString().trim().split("\t");
        final String[] tags = Arrays.copyOfRange(f, 11, f.length);
        Arrays.sort(tags);
        return String.join("\t", Arrays.copyOf(f, 11)) + "\t" + String.join("\t", tags);
    }

    @Test
    public void fixture_should_have_many_rans_embedded_reference_slices() throws IOException {
        int slices = 0;
        int ransEmbedded = 0;
        try (InputStream in = new FileInputStream(CRAM)) {
            final CramContainerIterator it = new CramContainerIterator(in);
            while (it.hasNext()) {
                final Container c = it.next();
                for (final Slice s : c.getSlices()) {
                    slices++;
                    if (s.getEmbeddedReferenceBlock() != null && s.getEmbeddedReferenceBlock().getCompressionMethod() == BlockCompressionMethod.RANS) {
                        ransEmbedded++;
                    }
                }
            }
        }
        Assert.assertTrue(slices + " slices, " + ransEmbedded + " rANS embedded", ransEmbedded > 5);
    }

    @Test
    public void decode_with_embedded_reference_should_match_htslib() throws IOException {
        final List<String> expected = new ArrayList<>();
        try (SamReader reader = SamReaderFactory.makeDefault().open(BAM)) {
            reader.forEach(r -> expected.add(fieldsAndSortedTags(r)));
        }
        Assert.assertEquals(330, expected.size());
        for (final int threads : new int[] { 0, 2 }) {
            final List<String> actual = new ArrayList<>();
            try (InputStream in = new FileInputStream(CRAM)) {
                final CRAMIterator it = new CRAMIterator(in, NO_REFERENCE, ValidationStringency.STRICT, threads, 2 * threads);
                it.forEachRemaining(r -> actual.add(fieldsAndSortedTags(r)));
                it.close();
            }
            Assert.assertEquals("threads " + threads, expected, actual);
        }
    }
}
