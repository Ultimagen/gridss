package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import org.junit.Test;

import static org.junit.Assert.*;

public class CoverageReadTest extends TestHelper {
	private final ReadPairConcordanceCalculator pairing = new FixedSizeReadPairConcordanceCalculator(0, 300);

	@Test
	public void isCounted_should_drop_unmapped_reads_and_duplicates_only_if_ignored() {
		SAMRecord dup = Read(0, 100, 10);
		dup.setDuplicateReadFlag(true);
		assertFalse(CoverageRead.isCounted(Unmapped(10), false));
		assertFalse(CoverageRead.isCounted(dup, true));
		assertTrue(CoverageRead.isCounted(dup, false));
		assertTrue(CoverageRead.isCounted(Read(0, 100, 10), true));
	}

	@Test
	public void create_should_record_the_mate_start_only_for_the_lower_read_of_a_non_overlapping_concordant_pair() {
		SAMRecord[] adjacent = RP(0, 100, 110, 10);
		assertEquals(110, CoverageRead.create(adjacent[0], pairing).pairMateStart);
		assertEquals(0, CoverageRead.create(adjacent[1], pairing).pairMateStart);
		SAMRecord[] overlapping = RP(0, 100, 109, 10);
		assertEquals(0, CoverageRead.create(overlapping[0], pairing).pairMateStart);
		SAMRecord[] discordant = RP(0, 100, 1000, 10);
		assertEquals(0, CoverageRead.create(discordant[0], pairing).pairMateStart);
	}
}
