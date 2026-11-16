package au.edu.wehi.idsv;

import com.google.common.collect.Iterators;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;

import java.io.Closeable;
import java.util.ArrayDeque;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.function.Supplier;

/**
 * One decode of a coordinate sorted input shared by the evidence iterator and the
 * reference coverage lookup of AnnotateVariants.
 *
 * Whichever consumer is ahead takes the next record from the file, under the tee's lock.
 * A record the evidence side takes is handed to the coverage side as a {@link CoverageRead},
 * captured before the evidence path mutates the record. A record the coverage side takes is
 * queued whole for the evidence side. Neither consumer ever waits for the other, so the tee
 * cannot deadlock the AnnotateVariants pipeline.
 *
 * Both queues are bounded. If either bound is exceeded, the coverage side detaches: it drains
 * the reads already queued for it, then decodes the file on its own from the first record it
 * has not seen. Its reads are then exactly those it would have had from the shared decode.
 */
public class CoverageReadTee {
	private static final Log log = Log.getInstance(CoverageReadTee.class);
	private static final int FIELDS = 4;
	private static final int CHUNK_RECORDS = 1024;
	private static final int EVIDENCE_BATCH = 64;
	private final String name;
	private final Supplier<SamReader> sharedReaderFactory;
	private final Supplier<SamReader> coverageReaderFactory;
	private final ReadPairConcordanceCalculator pairing;
	private final boolean ignoreDuplicates;
	private final int maxCoverageBuffered;
	private final int maxEvidencePending;
	// all state below is guarded by this
	private SamReader reader;
	private SAMRecordIterator raw;
	private boolean rawOpened;
	private long taken = 0;
	private long detachedAt = -1;
	private boolean evidenceClaimed;
	private boolean coverageClaimed;
	private boolean evidenceClosed;
	private boolean coverageClosed;
	private final ArrayDeque<SAMRecord> evidencePending = new ArrayDeque<>();
	private final ArrayDeque<Chunk> coverageChunks = new ArrayDeque<>();
	private Chunk tail;
	private int coverageBuffered = 0;
	private int maxCoverageBufferedSeen = 0;
	private int maxEvidencePendingSeen = 0;

	/**
	 * @param sharedReaderFactory opens the reader whose records both consumers see
	 * @param coverageReaderFactory opens the coverage side's own reader of the same file, used only after it detaches
	 * @param maxCoverageBuffered maximum reads queued for the coverage side
	 * @param maxEvidencePending maximum records queued for the evidence side
	 */
	public CoverageReadTee(String name, Supplier<SamReader> sharedReaderFactory, Supplier<SamReader> coverageReaderFactory,
			ReadPairConcordanceCalculator pairing, boolean ignoreDuplicates, int maxCoverageBuffered, int maxEvidencePending) {
		this.name = name;
		this.sharedReaderFactory = sharedReaderFactory;
		this.coverageReaderFactory = coverageReaderFactory;
		this.pairing = pairing;
		this.ignoreDuplicates = ignoreDuplicates;
		this.maxCoverageBuffered = maxCoverageBuffered;
		this.maxEvidencePending = maxEvidencePending;
	}

	/**
	 * Every record of the file, in file order, or null if already claimed.
	 */
	public synchronized CloseableIterator<SAMRecord> claimEvidenceIterator() {
		if (evidenceClaimed) return null;
		evidenceClaimed = true;
		return new EvidenceIterator();
	}

	/**
	 * Reads counted by {@link CoverageRead#isCounted(SAMRecord, boolean)}, in file order, or null if already claimed.
	 */
	public synchronized CloseableIterator<CoverageRead> claimCoverageIterator() {
		if (coverageClaimed) return null;
		coverageClaimed = true;
		return new CoverageIterator();
	}

	/**
	 * Raw record count at which the coverage side detached, -1 if it has not
	 */
	public synchronized long getDetachedAt() {
		return detachedAt;
	}

	private static final class Chunk {
		final int[] data = new int[CHUNK_RECORDS * FIELDS];
		int size = 0;
	}

	private SAMRecord takeRaw() {
		if (!rawOpened) {
			rawOpened = true;
			reader = sharedReaderFactory.get();
			raw = reader.iterator();
			raw.assertSorted(SortOrder.coordinate);
		}
		if (raw == null || !raw.hasNext()) return null;
		taken++;
		return raw.next();
	}

	private void offerToCoverage(SAMRecord r) {
		if (detachedAt >= 0 || coverageClosed) return;
		if (!CoverageRead.isCounted(r, ignoreDuplicates)) return;
		CoverageRead cr = CoverageRead.create(r, pairing);
		if (tail == null) {
			tail = new Chunk();
			coverageChunks.add(tail);
		}
		int offset = tail.size * FIELDS;
		tail.data[offset] = cr.referenceIndex;
		tail.data[offset + 1] = cr.alignmentStart;
		tail.data[offset + 2] = cr.alignmentEnd;
		tail.data[offset + 3] = cr.pairMateStart;
		tail.size++;
		if (tail.size == CHUNK_RECORDS) tail = null;
		coverageBuffered++;
		maxCoverageBufferedSeen = Math.max(maxCoverageBufferedSeen, coverageBuffered);
		if (coverageBuffered > maxCoverageBuffered) {
			detach(String.format("%d reads queued for coverage", coverageBuffered));
		}
	}

	private Chunk takeCoverageChunk() {
		Chunk c = coverageChunks.poll();
		if (c == tail) tail = null;
		if (c != null) coverageBuffered -= c.size;
		return c;
	}

	private void detach(String reason) {
		detachedAt = taken;
		log.warn(String.format("%s: coverage reads detached from the shared decode after %d records (%s). Coverage will decode the file again.", name, taken, reason));
	}

	private void closeSharedIfDone(String closing) {
		log.info(String.format("%s: %s closed after %d shared records; at most %d reads were queued for coverage, %d records for evidence; coverage detached at %d",
				name, closing, taken, maxCoverageBufferedSeen, maxEvidencePendingSeen, detachedAt));
		if (evidenceClosed && coverageClosed) {
			CloserUtil.close(raw);
			CloserUtil.close(reader);
			raw = null;
			reader = null;
			rawOpened = true;
		}
	}

	private class EvidenceIterator implements CloseableIterator<SAMRecord> {
		private final ArrayDeque<SAMRecord> batch = new ArrayDeque<>();
		@Override
		public boolean hasNext() {
			if (!batch.isEmpty()) return true;
			synchronized (CoverageReadTee.this) {
				if (evidenceClosed) return false;
				if (!evidencePending.isEmpty()) {
					batch.addAll(evidencePending);
					evidencePending.clear();
				} else {
					for (int i = 0; i < EVIDENCE_BATCH; i++) {
						SAMRecord r = takeRaw();
						if (r == null) break;
						offerToCoverage(r);
						batch.add(r);
					}
				}
			}
			return !batch.isEmpty();
		}
		@Override
		public SAMRecord next() {
			if (!hasNext()) throw new NoSuchElementException();
			return batch.poll();
		}
		@Override
		public void close() {
			synchronized (CoverageReadTee.this) {
				if (evidenceClosed) return;
				evidenceClosed = true;
				evidencePending.clear();
				closeSharedIfDone("evidence");
			}
		}
	}

	private class CoverageIterator implements CloseableIterator<CoverageRead> {
		private Chunk current;
		private int index;
		private CoverageRead pulled;
		private SamReader ownReader;
		private SAMRecordIterator ownRaw;
		private Iterator<CoverageRead> own;
		@Override
		public boolean hasNext() {
			if (pulled != null) return true;
			if (current != null && index < current.size) return true;
			if (own != null) return own.hasNext();
			current = null;
			long resumeAt;
			synchronized (CoverageReadTee.this) {
				if (coverageClosed) return false;
				Chunk c = takeCoverageChunk();
				if (c != null) {
					current = c;
					index = 0;
					return true;
				}
				while (detachedAt < 0) {
					if (evidencePending.size() >= maxEvidencePending) {
						detach(String.format("%d records queued for evidence", evidencePending.size()));
						break;
					}
					SAMRecord r = takeRaw();
					if (r == null) return false;
					if (!evidenceClosed) {
						evidencePending.add(r);
						maxEvidencePendingSeen = Math.max(maxEvidencePendingSeen, evidencePending.size());
					}
					if (CoverageRead.isCounted(r, ignoreDuplicates)) {
						pulled = CoverageRead.create(r, pairing);
						return true;
					}
				}
				resumeAt = detachedAt;
			}
			openOwn(resumeAt);
			return own.hasNext();
		}
		private void openOwn(long skip) {
			ownReader = coverageReaderFactory.get();
			ownRaw = ownReader.iterator();
			ownRaw.assertSorted(SortOrder.coordinate);
			for (long i = 0; i < skip && ownRaw.hasNext(); i++) {
				ownRaw.next();
			}
			own = Iterators.transform(Iterators.filter(ownRaw, r -> CoverageRead.isCounted(r, ignoreDuplicates)), r -> CoverageRead.create(r, pairing));
		}
		@Override
		public CoverageRead next() {
			if (!hasNext()) throw new NoSuchElementException();
			if (pulled != null) {
				CoverageRead r = pulled;
				pulled = null;
				return r;
			}
			if (current == null || index >= current.size) {
				return own.next();
			}
			int offset = index * FIELDS;
			index++;
			return new CoverageRead(current.data[offset], current.data[offset + 1], current.data[offset + 2], current.data[offset + 3]);
		}
		@Override
		public void close() {
			synchronized (CoverageReadTee.this) {
				if (coverageClosed) return;
				coverageClosed = true;
				coverageChunks.clear();
				tail = null;
				coverageBuffered = 0;
				closeSharedIfDone("coverage");
			}
			CloserUtil.close(ownRaw);
			CloserUtil.close(ownReader);
		}
	}
}
