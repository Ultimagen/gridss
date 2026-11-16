package htsjdk.samtools;

import htsjdk.samtools.cram.build.CramContainerIterator;
import htsjdk.samtools.cram.ref.CRAMReferenceSource;
import htsjdk.samtools.cram.ref.ReferenceContext;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.cram.structure.CRAMEncodingStrategy;
import htsjdk.samtools.cram.structure.Container;
import htsjdk.samtools.cram.structure.Slice;
import htsjdk.samtools.seekablestream.SeekableMemoryStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Random;
import java.util.function.Supplier;

public class CRAMIteratorParallelDecodeTest {
    private static final File REFERENCE = new File("src/test/resources/small.fa");
    private static final String DECODE_THREAD_NAME = "htsjdk-cram-decode";
    private static final int RECORD_COUNT = 1112;
    private static byte[] cram;
    private static SAMFileHeader header;

    @BeforeClass
    public static void writeCram() {
        header = new SAMFileHeader(SAMSequenceDictionaryExtractor.extractDictionary(REFERENCE.toPath()));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        header.addReadGroup(new SAMReadGroupRecord("rg"));
        cram = writeCram(createRecords(header));
    }

    /**
     * Small slices and containers so that the file has many containers, multi-ref slices for the contigs
     * with only a few reads, a container of placed unmapped reads and an unmapped-unplaced tail.
     */
    private static byte[] writeCram(final List<SAMRecord> records) {
        final CRAMEncodingStrategy strategy = new CRAMEncodingStrategy()
                .setMinimumSingleReferenceSliceSize(20)
                .setReadsPerSlice(50)
                .setSlicesPerContainer(2);
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        final CRAMFileWriter writer = new CRAMFileWriter(strategy, out, null, true, new ReferenceSource(REFERENCE), header, "test.cram");
        records.forEach(writer::addAlignment);
        writer.close();
        return out.toByteArray();
    }

    private static List<SAMRecord> createRecords(final SAMFileHeader header) {
        final Random random = new Random(0);
        // reads per contig: the 3-read contigs end up in multi-ref slices
        final int[] readsPerContig = { 400, 3, 3, 300, 3, 250, 3, 0 };
        final List<SAMRecord> records = new ArrayList<>();
        for (int contig = 0; contig < readsPerContig.length; contig++) {
            for (int i = 0; i < readsPerContig[contig]; i++) {
                final SAMRecord r = createRecord(header, random, "r" + records.size());
                r.setReferenceIndex(contig);
                r.setAlignmentStart(1 + random.nextInt(9000));
                r.setMappingQuality(random.nextInt(61));
                r.setCigarString(random.nextInt(4) == 0 ? "10S30M2D20M10I30M" : "100M");
                r.setReadPairedFlag(true);
                r.setFirstOfPairFlag(random.nextBoolean());
                r.setSecondOfPairFlag(!r.getFirstOfPairFlag());
                r.setReadNegativeStrandFlag(random.nextBoolean());
                r.setMateReferenceIndex(random.nextInt(header.getSequenceDictionary().size()));
                r.setMateAlignmentStart(1 + random.nextInt(9000));
                if (random.nextInt(10) == 0) {
                    // placed unmapped read, positioned at its mapped mate
                    r.setReadUnmappedFlag(true);
                    r.setCigarString("*");
                    r.setMappingQuality(0);
                }
                records.add(r);
            }
        }
        for (int i = 0; i < 150; i++) {
            final SAMRecord r = createRecord(header, random, "u" + i);
            r.setReadUnmappedFlag(true);
            records.add(r);
        }
        records.sort(Comparator.comparingInt((SAMRecord r) -> r.getReferenceIndex() < 0 ? Integer.MAX_VALUE : r.getReferenceIndex())
                .thenComparingInt(SAMRecord::getAlignmentStart));
        Assert.assertEquals(RECORD_COUNT, records.size());
        return records;
    }

    private static SAMRecord createRecord(final SAMFileHeader header, final Random random, final String name) {
        final SAMRecord r = new SAMRecord(header);
        r.setReadName(name);
        final byte[] bases = new byte[100];
        final byte[] quals = new byte[100];
        for (int i = 0; i < bases.length; i++) {
            bases[i] = (byte) "ACGT".charAt(random.nextInt(4));
            quals[i] = (byte) random.nextInt(41);
        }
        r.setReadBases(bases);
        r.setBaseQualities(quals);
        r.setAttribute("RG", "rg");
        r.setAttribute("t0", new String(quals).replaceAll("[^!-~]", "I"));
        r.setAttribute("tp", new byte[] { (byte) random.nextInt(3), 0, 1 });
        return r;
    }

    private static List<String> readAll(final CRAMIterator it) {
        final List<String> result = new ArrayList<>();
        try {
            while (it.hasNext()) {
                result.add(it.next().getSAMString());
            }
        } finally {
            it.close();
        }
        return result;
    }

    private static CRAMIterator streamIterator(final int threads, final int depth) {
        return new CRAMIterator(new ByteArrayInputStream(cram), new ReferenceSource(REFERENCE), ValidationStringency.STRICT, threads, depth);
    }

    private static CRAMIterator spanIterator(final QueryInterval[] intervals, final long[] coordinates, final int threads) {
        return new CRAMIterator(new SeekableMemoryStream(cram, "test.cram"), new ReferenceSource(REFERENCE), ValidationStringency.STRICT,
                intervals, coordinates, threads, 2 * threads);
    }

    private static List<Container> containers() {
        final List<Container> result = new ArrayList<>();
        final CramContainerIterator it = new CramContainerIterator(new ByteArrayInputStream(cram));
        while (it.hasNext()) {
            result.add(it.next());
        }
        return result;
    }

    private static long decodeThreadCount() {
        return Thread.getAllStackTraces().keySet().stream()
                .filter(t -> t.isAlive() && DECODE_THREAD_NAME.equals(t.getName()))
                .count();
    }

    private static void assertNoDecodeThreads() throws InterruptedException {
        for (int i = 0; i < 50 && decodeThreadCount() > 0; i++) {
            Thread.sleep(100);
        }
        Assert.assertEquals(0, decodeThreadCount());
    }

    @Test
    public void test_file_should_have_multiref_unmapped_and_many_containers() {
        final List<Container> containers = containers();
        Assert.assertTrue(String.valueOf(containers.size()), containers.size() > 8);
        final List<ReferenceContext> sliceContexts = new ArrayList<>();
        for (final Container c : containers) {
            for (final Slice s : c.getSlices()) {
                sliceContexts.add(s.getAlignmentContext().getReferenceContext());
            }
        }
        Assert.assertTrue(sliceContexts.stream().anyMatch(ReferenceContext::isMultiRef));
        Assert.assertTrue(sliceContexts.stream().anyMatch(ReferenceContext::isUnmappedUnplaced));
        Assert.assertTrue(sliceContexts.stream().anyMatch(ReferenceContext::isMappedSingleRef));
    }

    @Test
    public void parallel_decode_should_return_the_same_records_in_the_same_order() {
        final List<String> expected = readAll(streamIterator(0, 0));
        Assert.assertEquals(RECORD_COUNT, expected.size());
        for (final int[] threadsAndDepth : new int[][] { { 1, 1 }, { 1, 2 }, { 2, 2 }, { 4, 1 }, { 4, 8 }, { 8, 64 } }) {
            Assert.assertEquals(Arrays.toString(threadsAndDepth), expected, readAll(streamIterator(threadsAndDepth[0], threadsAndDepth[1])));
        }
    }

    @Test
    public void parallel_decode_should_match_serial_decode_for_query_spans() {
        final List<Container> containers = containers();
        final long[] coordinates = {
                containers.get(2).getContainerByteOffset() << 16, containers.get(6).getContainerByteOffset() << 16,
                containers.get(10).getContainerByteOffset() << 16, Long.MAX_VALUE };
        final QueryInterval[] intervals = QueryInterval.optimizeIntervals(new QueryInterval[] {
                new QueryInterval(0, 2000, 4000), new QueryInterval(3, 1, 5000), new QueryInterval(5, 8000, -1) });
        final List<String> expected = readAll(spanIterator(intervals, coordinates, 0));
        Assert.assertFalse(expected.isEmpty());
        Assert.assertEquals(expected, readAll(spanIterator(intervals, coordinates, 4)));
        Assert.assertEquals(readAll(spanIterator(null, coordinates, 0)), readAll(spanIterator(null, coordinates, 4)));
    }

    /**
     * CRAMFileReader.queryUnmapped() seeks the stream after constructing the iterator.
     */
    @Test
    public void parallel_decode_should_start_at_the_stream_position_of_the_first_read() throws IOException {
        final List<Container> containers = containers();
        final long tail = containers.get(containers.size() - 4).getContainerByteOffset();
        final List<List<String>> results = new ArrayList<>();
        for (final int threads : new int[] { 0, 4 }) {
            final SeekableStream stream = new SeekableMemoryStream(cram, "test.cram");
            final CRAMIterator it = new CRAMIterator(stream, new ReferenceSource(REFERENCE), ValidationStringency.STRICT, threads, 8);
            stream.seek(tail);
            boolean atAlignments;
            do {
                atAlignments = it.advanceToAlignmentInContainer(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX, SAMRecord.NO_ALIGNMENT_START);
            } while (!atAlignments && it.hasNext());
            results.add(readAll(it));
        }
        Assert.assertEquals(150, results.get(0).size());
        Assert.assertEquals(results.get(0), results.get(1));
    }

    @Test
    public void close_before_the_end_should_stop_decode_threads() throws InterruptedException {
        assertNoDecodeThreads();
        final CRAMIterator it = streamIterator(4, 8);
        for (int i = 0; i < 10; i++) {
            it.next();
        }
        Assert.assertTrue(decodeThreadCount() > 0);
        it.close();
        assertNoDecodeThreads();
    }

    @Test
    public void iterator_should_not_start_decode_threads_until_read() throws InterruptedException {
        assertNoDecodeThreads();
        final CRAMIterator it = streamIterator(4, 8);
        Assert.assertEquals(0, decodeThreadCount());
        it.close();
        assertNoDecodeThreads();
    }

    @Test
    public void reading_to_the_end_should_stop_decode_threads_without_close() throws InterruptedException {
        assertNoDecodeThreads();
        final CRAMIterator it = streamIterator(4, 8);
        int count = 0;
        while (it.hasNext()) {
            it.next();
            count++;
        }
        Assert.assertFalse(it.hasNext());
        Assert.assertEquals(RECORD_COUNT, count);
        assertNoDecodeThreads();
    }

    private static class DecodeFailure extends RuntimeException {
        DecodeFailure() {
            super("test decode failure");
        }
    }

    /**
     * Fails the decode of every slice on contig 3 (the 4th contig).
     */
    private static CRAMReferenceSource failingReference() {
        final ReferenceSource source = new ReferenceSource(REFERENCE);
        return new CRAMReferenceSource() {
            @Override
            public byte[] getReferenceBases(final SAMSequenceRecord sequenceRecord, final boolean tryNameVariants) {
                if (sequenceRecord.getSequenceIndex() == 3) throw new DecodeFailure();
                return source.getReferenceBases(sequenceRecord, tryNameVariants);
            }

            @Override
            public byte[] getReferenceBasesByRegion(final SAMSequenceRecord sequenceRecord, final int zeroBasedStart, final int requestedRegionLength) {
                if (sequenceRecord.getSequenceIndex() == 3) throw new DecodeFailure();
                return source.getReferenceBasesByRegion(sequenceRecord, zeroBasedStart, requestedRegionLength);
            }
        };
    }

    /**
     * @return records read before the exception, and the exception class
     */
    private static String readUntilException(final Supplier<CRAMIterator> iterator) {
        int count = 0;
        final CRAMIterator it = iterator.get();
        try {
            while (it.hasNext()) {
                it.next();
                count++;
            }
        } catch (final RuntimeException e) {
            return count + " " + e.getClass().getName();
        } finally {
            it.close();
        }
        throw new AssertionError("no exception after " + count + " records");
    }

    @Test
    public void decode_exception_should_be_rethrown_to_the_reader_after_the_preceding_records() throws InterruptedException {
        final String serial = readUntilException(() -> new CRAMIterator(new ByteArrayInputStream(cram), failingReference(), ValidationStringency.STRICT, 0, 0));
        Assert.assertTrue(serial, serial.endsWith(DecodeFailure.class.getName()));
        for (final int threads : new int[] { 1, 4 }) {
            Assert.assertEquals(serial, readUntilException(() -> new CRAMIterator(new ByteArrayInputStream(cram), failingReference(), ValidationStringency.STRICT, threads, 8)));
        }
        assertNoDecodeThreads();
    }

    @Test
    public void read_exception_should_be_rethrown_to_the_reader_after_the_preceding_records() throws InterruptedException {
        final List<Container> containers = containers();
        final byte[] truncated = Arrays.copyOf(cram, (int) containers.get(8).getContainerByteOffset() + 100);
        final Supplier<CRAMIterator> serial = () -> new CRAMIterator(new ByteArrayInputStream(truncated), new ReferenceSource(REFERENCE), ValidationStringency.STRICT, 0, 0);
        final Supplier<CRAMIterator> parallel = () -> new CRAMIterator(new ByteArrayInputStream(truncated), new ReferenceSource(REFERENCE), ValidationStringency.STRICT, 4, 8);
        Assert.assertEquals(readUntilException(serial), readUntilException(parallel));
        assertNoDecodeThreads();
    }
}
