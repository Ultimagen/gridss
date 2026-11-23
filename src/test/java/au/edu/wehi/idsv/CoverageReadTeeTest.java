package au.edu.wehi.idsv;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import org.junit.Test;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.function.Supplier;

import static org.junit.Assert.*;

public class CoverageReadTeeTest extends IntermediateFilesTest {
	private final ReadPairConcordanceCalculator pairing = new FixedSizeReadPairConcordanceCalculator(0, 300);
	private int opened = 0;
	private void createRandomInput(int pairs, long seed) {
		Random rng = new Random(seed);
		List<SAMRecord> reads = new ArrayList<>();
		for (int i = 0; i < pairs; i++) {
			int ref = rng.nextInt(3);
			int pos = 1 + rng.nextInt(5000);
			SAMRecord[] rp;
			switch (rng.nextInt(5)) {
				case 0:
					rp = DP(ref, pos, "10M5S", true, rng.nextInt(3), 1 + rng.nextInt(5000), "15M", false);
					break;
				case 1:
					rp = OEA(ref, pos, "15M", rng.nextBoolean());
					break;
				default:
					rp = RP(ref, pos, pos + rng.nextInt(400), 15);
					rp[0].setAttribute("MC", "15M");
					rp[1].setAttribute("MC", "15M");
					break;
			}
			for (SAMRecord r : rp) {
				r.setReadName(r.getReadName() + "#" + i);
				r.setDuplicateReadFlag(rng.nextInt(10) == 0);
				reads.add(r);
			}
		}
		createInput(reads);
	}
	private SamReader open() {
		opened++;
		return SamReaderFactory.makeDefault().open(input);
	}
	private CoverageReadTee tee(int maxCoverageBuffered, int maxEvidencePending) {
		Supplier<SamReader> factory = this::open;
		return new CoverageReadTee("test", factory, factory, pairing, true, maxCoverageBuffered, maxEvidencePending);
	}
	private List<String> expectedRecords() {
		try (SamReader reader = SamReaderFactory.makeDefault().open(input)) {
			List<String> result = new ArrayList<>();
			for (SAMRecord r : reader) result.add(r.getSAMString());
			return result;
		} catch (IOException e) {
			throw new UncheckedIOException(e);
		}
	}
	private List<String> expectedCoverage() {
		try (SamReader reader = SamReaderFactory.makeDefault().open(input)) {
			List<String> result = new ArrayList<>();
			for (SAMRecord r : reader) {
				if (CoverageRead.isCounted(r, true)) result.add(key(CoverageRead.create(r, pairing)));
			}
			return result;
		} catch (IOException e) {
			throw new UncheckedIOException(e);
		}
	}
	private static String key(CoverageRead r) {
		return r.referenceIndex + ":" + r.alignmentStart + "-" + r.alignmentEnd + ">" + r.pairMateStart;
	}
	private static List<String> drainEvidence(CloseableIterator<SAMRecord> it, boolean mutate) {
		List<String> result = new ArrayList<>();
		while (it.hasNext()) {
			SAMRecord r = it.next();
			result.add(r.getSAMString());
			if (mutate) {
				// what SAMEvidenceSource.transform can do to a record
				r.setReadUnmappedFlag(true);
				r.setMateUnmappedFlag(true);
				r.setAlignmentStart(1);
			}
		}
		it.close();
		return result;
	}
	private static List<String> drainCoverage(CloseableIterator<CoverageRead> it) {
		List<String> result = new ArrayList<>();
		while (it.hasNext()) result.add(key(it.next()));
		it.close();
		return result;
	}
	@Test
	public void coverage_should_see_records_as_decoded_when_evidence_runs_first_and_mutates_them() {
		createRandomInput(2000, 1);
		CoverageReadTee tee = tee(1 << 20, 1000);
		CloseableIterator<SAMRecord> eit = tee.claimEvidenceIterator();
		CloseableIterator<CoverageRead> cit = tee.claimCoverageIterator();
		assertEquals(expectedRecords(), drainEvidence(eit, true));
		assertEquals(expectedCoverage(), drainCoverage(cit));
		assertEquals(-1, tee.getDetachedAt());
		assertEquals(1, opened);
	}
	@Test
	public void evidence_should_get_all_records_when_coverage_runs_first() {
		createRandomInput(2000, 2);
		CoverageReadTee tee = tee(1 << 20, 1 << 20);
		CloseableIterator<SAMRecord> eit = tee.claimEvidenceIterator();
		CloseableIterator<CoverageRead> cit = tee.claimCoverageIterator();
		assertEquals(expectedCoverage(), drainCoverage(cit));
		assertEquals(expectedRecords(), drainEvidence(eit, true));
		assertEquals(-1, tee.getDetachedAt());
		assertEquals(1, opened);
	}
	@Test
	public void should_interleave_consumers_in_any_order() {
		createRandomInput(2000, 3);
		Random rng = new Random(3);
		CoverageReadTee tee = tee(1 << 20, 1 << 20);
		CloseableIterator<SAMRecord> eit = tee.claimEvidenceIterator();
		CloseableIterator<CoverageRead> cit = tee.claimCoverageIterator();
		List<String> evidence = new ArrayList<>();
		List<String> coverage = new ArrayList<>();
		while (eit.hasNext() || cit.hasNext()) {
			for (int i = rng.nextInt(200); i > 0 && eit.hasNext(); i--) {
				SAMRecord r = eit.next();
				evidence.add(r.getSAMString());
				r.setReadUnmappedFlag(true);
			}
			for (int i = rng.nextInt(200); i > 0 && cit.hasNext(); i--) {
				coverage.add(key(cit.next()));
			}
		}
		assertEquals(expectedRecords(), evidence);
		assertEquals(expectedCoverage(), coverage);
		assertEquals(1, opened);
	}
	@Test
	public void coverage_should_detach_and_decode_again_when_its_queue_is_full() {
		createRandomInput(2000, 4);
		CoverageReadTee tee = tee(100, 1000);
		CloseableIterator<SAMRecord> eit = tee.claimEvidenceIterator();
		CloseableIterator<CoverageRead> cit = tee.claimCoverageIterator();
		assertEquals(expectedRecords(), drainEvidence(eit, true));
		assertTrue(tee.getDetachedAt() > 0);
		assertEquals(expectedCoverage(), drainCoverage(cit));
		assertEquals(2, opened);
	}
	@Test
	public void coverage_should_detach_when_evidence_queue_is_full() {
		createRandomInput(2000, 5);
		CoverageReadTee tee = tee(1 << 20, 100);
		CloseableIterator<SAMRecord> eit = tee.claimEvidenceIterator();
		CloseableIterator<CoverageRead> cit = tee.claimCoverageIterator();
		assertEquals(expectedCoverage(), drainCoverage(cit));
		assertEquals(100, tee.getDetachedAt());
		assertEquals(expectedRecords(), drainEvidence(eit, true));
		assertEquals(2, opened);
	}
	@Test
	public void coverage_should_continue_after_evidence_closes() {
		createRandomInput(500, 6);
		CoverageReadTee tee = tee(1 << 20, 10);
		CloseableIterator<SAMRecord> eit = tee.claimEvidenceIterator();
		CloseableIterator<CoverageRead> cit = tee.claimCoverageIterator();
		eit.next();
		eit.close();
		assertEquals(expectedCoverage(), drainCoverage(cit));
		assertEquals(-1, tee.getDetachedAt());
	}
	@Test
	public void should_claim_each_consumer_once() {
		createRandomInput(10, 7);
		CoverageReadTee tee = tee(1 << 20, 1000);
		assertNotNull(tee.claimEvidenceIterator());
		assertNull(tee.claimEvidenceIterator());
		assertNotNull(tee.claimCoverageIterator());
		assertNull(tee.claimCoverageIterator());
	}
	@Test
	public void should_match_independent_decodes_on_concurrent_consumers() throws Exception {
		createRandomInput(20000, 8);
		List<String> expectedRecords = expectedRecords();
		List<String> expectedCoverage = expectedCoverage();
		ExecutorService pool = Executors.newFixedThreadPool(2);
		try {
			for (int[] limits : new int[][] { { 1 << 20, 1 << 20 }, { 500, 1 << 20 }, { 1 << 20, 50 } }) {
				for (int round = 0; round < 3; round++) {
					CoverageReadTee tee = tee(limits[0], limits[1]);
					CloseableIterator<SAMRecord> eit = tee.claimEvidenceIterator();
					CloseableIterator<CoverageRead> cit = tee.claimCoverageIterator();
					Future<List<String>> evidence = pool.submit(() -> drainEvidence(eit, true));
					Future<List<String>> coverage = pool.submit(() -> drainCoverage(cit));
					assertEquals(expectedRecords, evidence.get());
					assertEquals(expectedCoverage, coverage.get());
				}
			}
		} finally {
			pool.shutdownNow();
		}
	}
	@Test
	public void lookup_from_coverage_reads_should_match_lookup_from_records() {
		createRandomInput(3000, 9);
		List<SAMRecord> records;
		try (SamReader reader = SamReaderFactory.makeDefault().open(input)) {
			records = Lists.newArrayList(reader.iterator());
		} catch (IOException e) {
			throw new UncheckedIOException(e);
		}
		CoverageReadTee tee = tee(1 << 20, 1000);
		CloseableIterator<SAMRecord> eit = tee.claimEvidenceIterator();
		drainEvidence(eit, true);
		ReferenceCoverageLookup fromRecords = new SequentialReferenceCoverageLookup(records.iterator(), IDSV(records), pairing, 1000, 0, true);
		ReferenceCoverageLookup fromTee = SequentialReferenceCoverageLookup.fromCoverageReads(tee.claimCoverageIterator(), IDSV(records), pairing, 1000, 0);
		for (int ref = 0; ref < 3; ref++) {
			for (int pos = 1; pos < 5500; pos += 7) {
				assertEquals(fromRecords.readsSupportingNoBreakendAfter(ref, pos), fromTee.readsSupportingNoBreakendAfter(ref, pos));
				assertEquals(fromRecords.readPairsSupportingNoBreakendAfter(ref, pos), fromTee.readPairsSupportingNoBreakendAfter(ref, pos));
			}
		}
	}
	@Test
	public void lookup_should_match_when_reads_outside_countable_windows_are_dropped() {
		createRandomInput(3000, 10);
		List<SAMRecord> records;
		try (SamReader reader = SamReaderFactory.makeDefault().open(input)) {
			records = Lists.newArrayList(reader.iterator());
		} catch (IOException e) {
			throw new UncheckedIOException(e);
		}
		int windowSize = 50;
		int maxEvidenceWindow = SequentialReferenceCoverageLookup.maxEvidenceWindow(IDSV(records), pairing);
		List<int[]> queries = Lists.newArrayList(new int[] { 0, 800, 820 }, new int[] { 0, 3000, 3005 }, new int[] { 0, 3010, 3012 }, new int[] { 2, 4000, 4100 });
		CoverageReadTee counting = tee(1 << 20, 1000);
		counting.setCountableWindows(new CoverageReadWindows(queries, windowSize, maxEvidenceWindow, 3));
		drainEvidence(counting.claimEvidenceIterator(), true);
		int queued = drainCoverage(counting.claimCoverageIterator()).size();
		assertTrue(queued > 0);
		assertTrue(queued < expectedCoverage().size() / 2);

		CoverageReadTee tee = tee(1 << 20, 1000);
		tee.setCountableWindows(new CoverageReadWindows(queries, windowSize, maxEvidenceWindow, 3));
		drainEvidence(tee.claimEvidenceIterator(), true);
		ReferenceCoverageLookup fromRecords = new SequentialReferenceCoverageLookup(records.iterator(), IDSV(records), pairing, windowSize, 0, true);
		ReferenceCoverageLookup fromTee = SequentialReferenceCoverageLookup.fromCoverageReads(tee.claimCoverageIterator(), IDSV(records), pairing, windowSize, 0);
		for (int[] q : queries) {
			for (int pos = q[1]; pos <= q[2]; pos++) {
				assertEquals(fromRecords.readsSupportingNoBreakendAfter(q[0], pos), fromTee.readsSupportingNoBreakendAfter(q[0], pos));
				assertEquals(fromRecords.readPairsSupportingNoBreakendAfter(q[0], pos), fromTee.readPairsSupportingNoBreakendAfter(q[0], pos));
			}
		}
	}
	@Test
	public void countable_windows_should_extend_before_each_query_by_both_windows() {
		CoverageReadWindows w = new CoverageReadWindows(Lists.newArrayList(new int[] { 0, 100, 110 }, new int[] { 0, 105, 200 }, new int[] { 1, 50, 50 }), 10, 5, 3);
		assertFalse(w.isCountable(0, 83));
		assertTrue(w.isCountable(0, 84));
		assertTrue(w.isCountable(0, 200));
		assertFalse(w.isCountable(0, 201));
		assertTrue(w.isCountable(1, 34));
		assertTrue(w.isCountable(1, 50));
		assertFalse(w.isCountable(1, 51));
		assertFalse(w.isCountable(2, 1));
	}
}
