/*******************************************************************************
 * Copyright 2013-2016 EMBL-EBI
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License countingInputStream distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ******************************************************************************/
package htsjdk.samtools;

import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.cram.build.CRAMReferenceRegion;
import htsjdk.samtools.cram.build.CramContainerIterator;
import htsjdk.samtools.cram.build.CramSpanContainerIterator;
import htsjdk.samtools.cram.io.CountingInputStream;
import htsjdk.samtools.cram.ref.CRAMReferenceSource;
import htsjdk.samtools.cram.structure.*;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.AsyncReadTaskRunner;
import htsjdk.samtools.util.Tuple;

import java.io.Closeable;
import java.io.IOException;
import java.io.InputStream;
import java.util.*;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import htsjdk.samtools.util.RuntimeIOException;

/**
 * GRIDSS fork of the htsjdk 3.0.1 CRAMIterator.
 *
 * With samjdk.cram_decode_threads=N (N > 0), containers are read on one thread and decoded on N threads
 * owned by this iterator, at most samjdk.cram_decode_depth containers ahead of the reader, and returned in
 * file order. Decode is the stock Container.getSAMRecords(), with one reference region and compressor
 * cache per decode thread. With N = 0 (the default) the iterator behaves exactly as in htsjdk.
 */
public class CRAMIterator implements SAMRecordIterator, Closeable {
    private final CountingInputStream countingInputStream;
    private final CramContainerIterator containerIterator;
    private final CramHeader cramHeader;
    private final SAMFileHeader samFileHeader;
    private final CRAMReferenceRegion cramReferenceRegion;
    private final QueryInterval[] queryIntervals;
    private final CRAMReferenceSource referenceSource;
    private final int decodeThreads;
    private final int decodeDepth;
    private ThreadPoolExecutor decodePool;
    private ContainerDecoder containerDecoder;

    private volatile ValidationStringency validationStringency;
    private List<SAMRecord> samRecords;
    private Container container;
    private SamReader mReader;
    private final long firstContainerOffset;

    // Keep a cache of re-usable compressor instances to reduce the need to repeatedly reallocate
    // large numbers of small temporary objects, especially for the RANS compressor, which
    // allocates ~256k small objects every time its instantiated.
    private final CompressorCache compressorCache = new CompressorCache();

    /**
     * `samRecordIndex` only used when validation is not `SILENT`
     * (for identification by the validator which records are invalid)
     */
    private long samRecordIndex;
    private Iterator<SAMRecord> samRecordIterator = Collections.EMPTY_LIST.iterator();

    public CRAMIterator(final InputStream inputStream,
                        final CRAMReferenceSource referenceSource,
                        final ValidationStringency validationStringency) {
        this(inputStream, referenceSource, validationStringency, Defaults.CRAM_DECODE_THREADS, Defaults.CRAM_DECODE_DEPTH);
    }

    CRAMIterator(final InputStream inputStream,
                 final CRAMReferenceSource referenceSource,
                 final ValidationStringency validationStringency,
                 final int decodeThreads,
                 final int decodeDepth) {
        this.countingInputStream = new CountingInputStream(inputStream);
        this.containerIterator = new CramContainerIterator(this.countingInputStream);

        this.validationStringency = validationStringency;
        samFileHeader = containerIterator.getSamFileHeader();
        cramReferenceRegion = new CRAMReferenceRegion(referenceSource, samFileHeader.getSequenceDictionary());
        cramHeader = containerIterator.getCramHeader();
        firstContainerOffset = this.countingInputStream.getCount();
        samRecords = new ArrayList<>(new CRAMEncodingStrategy().getReadsPerSlice());
        this.queryIntervals = null;
        this.referenceSource = referenceSource;
        this.decodeThreads = decodeThreads;
        this.decodeDepth = Math.max(1, decodeDepth);
    }

    public CRAMIterator(final SeekableStream seekableStream,
                        final CRAMReferenceSource referenceSource,
                        final ValidationStringency validationStringency,
                        final QueryInterval[] queryIntervals,
                        final long[] coordinates) {
        this(seekableStream, referenceSource, validationStringency, queryIntervals, coordinates,
                Defaults.CRAM_DECODE_THREADS, Defaults.CRAM_DECODE_DEPTH);
    }

    CRAMIterator(final SeekableStream seekableStream,
                 final CRAMReferenceSource referenceSource,
                 final ValidationStringency validationStringency,
                 final QueryInterval[] queryIntervals,
                 final long[] coordinates,
                 final int decodeThreads,
                 final int decodeDepth) {
        this.countingInputStream = new CountingInputStream(seekableStream);
        this.containerIterator = CramSpanContainerIterator.fromFileSpan(seekableStream, coordinates);

        this.validationStringency = validationStringency;
        samFileHeader = containerIterator.getSamFileHeader();
        cramReferenceRegion = new CRAMReferenceRegion(referenceSource, samFileHeader.getSequenceDictionary());
        cramHeader = containerIterator.getCramHeader();
        firstContainerOffset = this.countingInputStream.getCount();
        samRecords = new ArrayList<>(new CRAMEncodingStrategy().getReadsPerSlice());
        this.queryIntervals = queryIntervals;
        this.referenceSource = referenceSource;
        this.decodeThreads = decodeThreads;
        this.decodeDepth = Math.max(1, decodeDepth);
    }

    private BAMIteratorFilter.FilteringIteratorState nextContainer() {
        if (decodeThreads > 0) {
            return nextDecodedContainer();
        }
        if (containerIterator != null) {
            if (!containerIterator.hasNext()) {
                samRecords.clear();
                return BAMIteratorFilter.FilteringIteratorState.STOP_ITERATION;
            }
            container = containerIterator.next();
            if (container.isEOF()) {
                samRecords.clear();
                return BAMIteratorFilter.FilteringIteratorState.STOP_ITERATION;
            }
        } else {
            final long containerByteOffset = countingInputStream.getCount();
            container = new Container(cramHeader.getCRAMVersion(), countingInputStream, containerByteOffset);
            if (container.isEOF()) {
                samRecords.clear();
                return BAMIteratorFilter.FilteringIteratorState.STOP_ITERATION;
            }
        }

        if (containerMatchesQuery(container)) {
            samRecords = container.getSAMRecords(
                    validationStringency,
                    cramReferenceRegion,
                    compressorCache,
                    getSAMFileHeader());
            samRecordIterator = samRecords.iterator();
            return BAMIteratorFilter.FilteringIteratorState.MATCHES_FILTER;
        } else {
            return BAMIteratorFilter.FilteringIteratorState.CONTINUE_ITERATION;
        }
    }

    private BAMIteratorFilter.FilteringIteratorState nextDecodedContainer() {
        // Started on first use, not in the constructor: CRAMFileReader.queryUnmapped() seeks the stream after
        // constructing the iterator, and CRAMFileReader keeps a header-only iterator that is never read.
        if (containerDecoder == null) {
            decodePool = new ThreadPoolExecutor(decodeThreads, decodeThreads, 60, TimeUnit.SECONDS,
                    new LinkedBlockingQueue<>(), r -> {
                        final Thread t = new Thread(r, "htsjdk-cram-decode");
                        t.setDaemon(true);
                        return t;
                    });
            // idle threads of an iterator that is neither exhausted nor closed must not live forever
            decodePool.allowCoreThreadTimeOut(true);
            containerDecoder = new ContainerDecoder(decodePool);
        }
        final List<SAMRecord> records;
        try {
            records = containerDecoder.nextRecord();
        } catch (final IOException e) {
            throw new RuntimeIOException(e);
        }
        if (records == null) {
            stopDecoding();
            samRecords.clear();
            return BAMIteratorFilter.FilteringIteratorState.STOP_ITERATION;
        }
        samRecords = records;
        samRecordIterator = samRecords.iterator();
        return BAMIteratorFilter.FilteringIteratorState.MATCHES_FILTER;
    }

    private void stopDecoding() {
        if (containerDecoder != null) {
            containerDecoder.disableAsyncProcessing();
            // waits for in-flight tasks, so the pool is idle when it is shut down
            containerDecoder.flushAsyncProcessing();
            decodePool.shutdown();
        }
    }

    /**
     * Reads containers serially (blocking I/O pool) and decodes each one as a separate batch on this iterator's
     * decode pool. AsyncReadTaskRunner returns batches in read order and bounds them to decodeDepth.
     */
    private final class ContainerDecoder extends AsyncReadTaskRunner<List<SAMRecord>, Container> {
        private final CRAMReferenceSource lockedReferenceSource = new LockedReferenceSource(referenceSource);
        private final ThreadLocal<CRAMReferenceRegion> referenceRegion = ThreadLocal.withInitial(() ->
                new CRAMReferenceRegion(lockedReferenceSource, samFileHeader.getSequenceDictionary()));
        private final ThreadLocal<CompressorCache> compressorCache = ThreadLocal.withInitial(CompressorCache::new);

        ContainerDecoder(final ThreadPoolExecutor pool) {
            super(1, decodeDepth, pool);
        }

        @Override
        public Tuple<Container, Long> performReadAhead(final long bufferBudget) {
            while (containerIterator.hasNext()) {
                final Container next = containerIterator.next();
                if (next.isEOF()) {
                    return null;
                }
                if (containerMatchesQuery(next)) {
                    return new Tuple<>(next, 1L);
                }
            }
            return null;
        }

        @Override
        public List<SAMRecord> transform(final Container record) {
            return record.getSAMRecords(validationStringency, referenceRegion.get(), compressorCache.get(), samFileHeader);
        }
    }

    /**
     * ReferenceSource is not thread-safe. Locks on the wrapped source, which iterators of the same reader share.
     */
    private static final class LockedReferenceSource implements CRAMReferenceSource {
        private final CRAMReferenceSource source;

        LockedReferenceSource(final CRAMReferenceSource source) {
            this.source = source;
        }

        @Override
        public byte[] getReferenceBases(final SAMSequenceRecord sequenceRecord, final boolean tryNameVariants) {
            synchronized (source) {
                return source.getReferenceBases(sequenceRecord, tryNameVariants);
            }
        }

        @Override
        public byte[] getReferenceBasesByRegion(final SAMSequenceRecord sequenceRecord, final int zeroBasedStart, final int requestedRegionLength) {
            synchronized (source) {
                return source.getReferenceBasesByRegion(sequenceRecord, zeroBasedStart, requestedRegionLength);
            }
        }
    }

    private boolean containerMatchesQuery(final Container container) {
        if (queryIntervals == null) {
            return true;
        } else {
            // binary search our query intervals to see if the alignment span of this container
            // overlaps any query - it doesn't matter which one, we only care whether or not there is a match
            final AlignmentContext alignmentContext = container.getAlignmentContext();
            return (!alignmentContext.getReferenceContext().isMappedSingleRef() ||
                Arrays.binarySearch(
                    queryIntervals,
                    new QueryInterval(
                            alignmentContext.getReferenceContext().getReferenceContextID(),
                            alignmentContext.getAlignmentStart(),
                            alignmentContext.getAlignmentStart() + alignmentContext.getAlignmentSpan() - 1
                    ),
                    overlapsContainerSpan) >= 0);
        }
    }

    //TODO: this should filter at the slice level!
    //we don't actually care which QueryInterval overlaps with the container; we just want to know if there is one...
    private final static Comparator<QueryInterval> overlapsContainerSpan = (queryInterval, containerInterval) -> {
        int comp = queryInterval.referenceIndex - containerInterval.referenceIndex;
        if (comp != 0) {
            return comp;
        }
        if (queryInterval.end <= 0) {
            // our query interval specifies a symbolic end, so call it a match if the container span
            // overlaps the start of the queryInterval
            return containerInterval.end <= queryInterval.start ?
                    -1 :
                    0;
        } else if (containerInterval.overlaps(queryInterval)) {
            return 0; // there is overlap so call it a match
        }
        return queryInterval.compareTo(containerInterval);
    };

    /**
     * Skip cached records until given alignment start position.
     *
     * @param refIndex reference sequence index
     * @param pos      alignment start to skip to
     */
    //TODO: this should first select the correct slice so we don't decode all slices unnecessarily
    public boolean advanceToAlignmentInContainer(final int refIndex, final int pos) {
        if (!hasNext()) return false;
        int i = 0;
        for (final SAMRecord record : samRecords) {
            if (refIndex != SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && record.getReferenceIndex() != refIndex) {
                continue;
            }

            if (pos <= 0) {
                if (record.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START) {
                    samRecordIterator = samRecords.listIterator(i);
                    return true;
                }
            } else {
                if (record.getAlignmentStart() >= pos) {
                    samRecordIterator = samRecords.listIterator(i);
                    return true;
                }
            }
            i++;
        }
        samRecordIterator = Collections.EMPTY_LIST.iterator();
        return false;
    }

    @Override
    public boolean hasNext() {
        if (container != null && container.isEOF()) {
            return false;
        }

        if (!samRecordIterator.hasNext()) {
            BAMIteratorFilter.FilteringIteratorState nextContainerPasses =
                    BAMIteratorFilter.FilteringIteratorState.CONTINUE_ITERATION;
            while (nextContainerPasses == BAMIteratorFilter.FilteringIteratorState.CONTINUE_ITERATION){
                nextContainerPasses = nextContainer();
            }
            return nextContainerPasses == BAMIteratorFilter.FilteringIteratorState.MATCHES_FILTER;
        }

        return !samRecords.isEmpty();
    }

    @Override
    public SAMRecord next() {
        if (hasNext()) {
            SAMRecord samRecord = samRecordIterator.next();
            if (validationStringency != ValidationStringency.SILENT) {
                SAMUtils.processValidationErrors(samRecord.isValid(), samRecordIndex++, validationStringency);
            }
            return samRecord;
        } else {
            throw new NoSuchElementException();
        }
    }

    @Override
    public void remove() {
        throw new RuntimeException("Removal of records not implemented.");
    }

    @Override
    public void close() {
        stopDecoding();
        samRecords.clear();
        try {
            if (countingInputStream != null) {
                countingInputStream.close();
            }
        } catch (final RuntimeIOException e) { }
    }

    public long getFirstContainerOffset() {
        return firstContainerOffset;
    }

    @Override
    public SAMRecordIterator assertSorted(final SortOrder sortOrder) {
        return SamReader.AssertingIterator.of(this).assertSorted(sortOrder);
    }

    public CramHeader getCramHeader() {
        return cramHeader;
    }

    public ValidationStringency getValidationStringency() {
        return validationStringency;
    }

    public void setValidationStringency(final ValidationStringency validationStringency) {
        this.validationStringency = validationStringency;
    }

    public SamReader getFileSource() {
        return mReader;
    }

    public void setFileSource(final SamReader mReader) {
        this.mReader = mReader;
    }

    public SAMFileHeader getSAMFileHeader() {
        return samFileHeader;
    }

}
