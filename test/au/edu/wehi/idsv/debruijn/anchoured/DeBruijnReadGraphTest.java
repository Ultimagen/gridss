package au.edu.wehi.idsv.debruijn.anchoured;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import htsjdk.samtools.SAMRecord;

import org.junit.Test;

import au.edu.wehi.idsv.BreakendDirection;
import au.edu.wehi.idsv.debruijn.anchored.DeBruijnReadGraph;
import au.edu.wehi.idsv.sam.AnomolousReadAssembly;
import au.edu.wehi.idsv.TestHelper;

public class DeBruijnReadGraphTest extends TestHelper {
	private static SAMRecord R(String read) {
		return R(null, read, null, false, true);
	}
	private static SAMRecord R(String readName, String read) {
		return R(readName, read, null, false, true);
	}
	private static SAMRecord R(String readName, String read, byte[] qual, boolean mappedNegativeStrand, boolean mateNegativeStrand) {
		SAMRecord record = new SAMRecord(getHeader());
		if (qual == null) {
			qual = new byte[read.length()];
			for (int i = 0; i < qual.length; i++) qual[i] = 1;
		}
		record.setReadBases(B(read));
		record.setBaseQualities(qual);
		record.setReadPairedFlag(true);
		record.setReadNegativeStrandFlag(mappedNegativeStrand);
		record.setReadUnmappedFlag(false);
		record.setMateUnmappedFlag(false);
		record.setMateNegativeStrandFlag(mateNegativeStrand);
		if (readName == null) {
			readName = String.format("%s-%s-%s%s", read, qual, mappedNegativeStrand, mateNegativeStrand);
		}
		record.setReadName(readName);
		return record;
	}
	@Test
	public void should_assemble_single_read() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(4, BreakendDirection.Forward);
		ass.addRead(R("AAAACGTC"), true);
		SAMRecord result = ass.assembleVariant();
		assertEquals("AAAACGTC", S(result.getReadBases()));
	}
	@Test
	public void should_assemble_positive_strand_consensus() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(4, BreakendDirection.Forward);
		ass.addRead(R(null, "AAAACGTC", null, true, true), true);
		SAMRecord result = ass.assembleVariant();
		assertEquals("AAAACGTC", S(result.getReadBases()));
	}
	@Test
	public void should_assemble_unanchored_reads() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(3, BreakendDirection.Forward);
		ass.addRead(R("CTAAA"), true);
		ass.addRead(R(null, "AAAGT", null, false, true), false);
		SAMRecord result = ass.assembleVariant();
		assertEquals("CTAAAGT", S(result.getReadBases()));
		assertEquals("5M2S", result.getCigarString());
	}
	@Test(expected = RuntimeException.class)  
	public void unanchored_reads_should_require_mapped_mate() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(3, BreakendDirection.Forward);
		ass.addRead(R("CTAAA"), true);
		SAMRecord unanchored = R(null, "AAAGT", null, false, true);
		unanchored.setMateUnmappedFlag(true);
		ass.addRead(unanchored, false);
		SAMRecord result = ass.assembleVariant();
		assertEquals("CTAAAGT", S(result.getReadBases()));
		assertEquals("5M2S", result.getCigarString());
	}
	@Test
	public void should_assemble_unanchored_reads_in_FR_orientation() {
		assertExpected(BreakendDirection.Forward, "CTAAAGT", "5M2S", "CTAAA", "AAAGT", true, false);
		assertExpected(BreakendDirection.Forward, "CTAAAGT", "5M2S", "CTAAA", "AAAGT", false, true);
		assertExpected(BreakendDirection.Forward, "CTAAAGT", "5M2S", "CTAAA", "ACTTT", true, true);
		assertExpected(BreakendDirection.Forward, "CTAAAGT", "5M2S", "CTAAA", "ACTTT", false, false);
		
		assertExpected(BreakendDirection.Backward, "GTAAACT", "2S5M", "AAACT", "GTAAA", true, false);
		assertExpected(BreakendDirection.Backward, "GTAAACT", "2S5M", "AAACT", "GTAAA", false, true);
		assertExpected(BreakendDirection.Backward, "GTAAACT", "2S5M", "AAACT", "TTTAC", true, true);
		assertExpected(BreakendDirection.Backward, "GTAAACT", "2S5M", "AAACT", "TTTAC", false, false);
	}
	private void assertExpected(BreakendDirection direction, String expectedSeq, String expectedCigar , String anchorSeq, String unanchorSeq, boolean mappedNegativeStrand, boolean mateNegativeStrand) {
		// Assembly should not depend on whether the read is mapped or not 
		DeBruijnReadGraph ass = new DeBruijnReadGraph(3, direction);
		ass.addRead(R(anchorSeq), true);
		SAMRecord unanchored = R(null, unanchorSeq, null, mappedNegativeStrand, mateNegativeStrand);
		unanchored.setReadUnmappedFlag(true);
		ass.addRead(unanchored, false);
		SAMRecord result = ass.assembleVariant();
		assertEquals(expectedSeq, S(result.getReadBases()));
		assertEquals(expectedCigar, result.getCigarString());
		
		ass = new DeBruijnReadGraph(3, direction);
		ass.addRead(R(anchorSeq), true);
		unanchored = R(null, unanchorSeq, null, mappedNegativeStrand, mateNegativeStrand);
		unanchored.setReadUnmappedFlag(false);
		ass.addRead(unanchored, false);
		result = ass.assembleVariant();
		assertEquals(expectedSeq, S(result.getReadBases()));
		assertEquals(expectedCigar, result.getCigarString());
	}
	@Test
	public void should_assemble_soft_clipped_read() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(4, BreakendDirection.Forward);
		SAMRecord sc = R("AAAACGTC");
		sc.setCigarString("4M4S");
		ass.addRead(sc, true);
		SAMRecord result = ass.assembleVariant();
		assertEquals("AAAACGTC", S(result.getReadBases()));
		assertEquals("4M4S", result.getCigarString());
	}
	@Test
	public void should_greedy_traverse_highest_weight_path() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(2, BreakendDirection.Forward);
		SAMRecord sc = R(null, "ACGTACTGAG", new byte[] { 1,2,3,4,5,6,7,8,9,10}, false, true);
		sc.setCigarString("4M6S");
		ass.addRead(sc, true);
		SAMRecord result = ass.assembleVariant();
		assertEquals("TACTGAGT", S(result.getReadBases()));
	}
	@Test
	public void should_assemble_backward_breakpoint() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(4, BreakendDirection.Backward);
		SAMRecord sc = R("AAAACGTC");
		sc.setCigarString("4S4M");
		ass.addRead(sc, true);
		SAMRecord result = ass.assembleVariant();
		assertEquals("AAAACGTC", S(result.getReadBases()));
		assertEquals("4S4M", result.getCigarString());
	}
	@Test
	public void remove_should_exclude_from_graph() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(4, BreakendDirection.Forward);
		SAMRecord sc = R("AAAACGTC");
		sc.setCigarString("4M4S");
		ass.addRead(sc, true);
		ass.addRead(R("r1", "AAAATTTT"), false);
		ass.addRead(R("r2", "AAAATTTT"), false);
		assertEquals("AAAATTTT", S(ass.assembleVariant().getReadBases()));
		ass.removeRead(R("r1", "AAAATTTT"), false);
		ass.removeRead(R("r2", "AAAATTTT"), false);
		SAMRecord result = ass.assembleVariant();
		assertEquals("AAAACGTC", S(result.getReadBases()));
	}
	@Test
	public void should_use_offset_kmer_if_softclip_longer_than_k() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(4, BreakendDirection.Forward);
		SAMRecord sc = R("AAAACGTC");
		sc.setCigarString("2M6S");
		ass.addRead(sc, true);
		SAMRecord result = ass.assembleVariant();
		assertEquals("AAAACGTC", S(result.getReadBases()));
		assertEquals("2M6S", result.getCigarString());
	}
	@Test
	public void assembly_base_quality_should_be_sum_of_min_read_qualities() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(3, BreakendDirection.Forward);
		ass.addRead(R(null, "ACGTA", new byte[] { 1,2,3,4,5 }, false, true), true);
		ass.addRead(R(null, "ACGTA", new byte[] { 3,4,5,6,7 }, false, true), true);
		SAMRecord result = ass.assembleVariant();
		// pad out read qualities
		assertArrayEquals(new byte[] { 4,6,8,8,8 }, result.getBaseQualities());
	}
	@Test
	public void assembly_base_quality_should_pad_to_match_read_length() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(3, BreakendDirection.Forward);
		ass.addRead(R(null, "ACGTA", new byte[] { 1,2,3,4,5 }, false, true), true);
		ass.addRead(R(null, "ACGTA", new byte[] { 3,4,5,6,7 }, false, true), true);
		SAMRecord result = ass.assembleVariant();
		// pad out read qualities
		assertEquals(result.getReadBases().length, result.getBaseQualities().length);
	}
	@Test
	public void read_count_should_be_number_of_reads__with_at_least_one_kmer_on_path() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(3, BreakendDirection.Forward);
		ass.addRead(R(null, "ACGTT", new byte[] { 1,2,3,4,5 }, false, true), true);
		ass.addRead(R(null, "ACGTA", new byte[] { 3,4,5,6,7 }, false, true), true);
		ass.addRead(R(null, "AAAAA", new byte[] { 1,1,1,1,1 }, false, true), true);
		ass.addRead(R(null, "ACGAA", new byte[] { 1,1,1,1,1 }, false, true), true);
		AnomolousReadAssembly result = ass.assembleVariant();
		
		assertEquals(3, (int)result.getReadCount());
	}
	@Test
	public void read_base_count_should_be_number_of_read_bases_on_returned_path() {
		DeBruijnReadGraph ass = new DeBruijnReadGraph(3, BreakendDirection.Forward);
		ass.addRead(R(null, "ACGTT", new byte[] { 1,2,3,4,5 }, false, true), true); // 4
		ass.addRead(R(null, "ACGTA", new byte[] { 3,4,5,6,7 }, false, true), true); // 5
		ass.addRead(R(null, "AAAAA", new byte[] { 1,1,1,1,1 }, false, true), true); // 0 since no kmer on path
		ass.addRead(R(null, "ACGAA", new byte[] { 1,1,1,1,1 }, false, true), true); // 3
		ass.addRead(R(null, "TTGTA", new byte[] { 1,1,1,1,1 }, false, true), false); // 3
		ass.addRead(R(null, "GTACG", new byte[] { 1,1,1,1,1 }, false, true), false); // 4+3 since kmers are disconnected 
		AnomolousReadAssembly result = ass.assembleVariant();
		
		assertEquals("test assumes this contruction - did I get the test case wrong?", "TACGTA", S(result.getReadBases()));
		assertEquals(4 + 5 + 3 + 3 + 4+3, (int)result.getReadBaseCount());
	}
}