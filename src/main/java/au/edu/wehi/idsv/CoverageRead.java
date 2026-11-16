package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;

/**
 * The fields of an aligned read that reference coverage counting uses.
 *
 * Taken from the SAMRecord as it was decoded, so the evidence path is free to
 * mutate the record afterwards.
 */
public final class CoverageRead {
	public final int referenceIndex;
	public final int alignmentStart;
	public final int alignmentEnd;
	/**
	 * Mate alignment start if this read is the lower mapped read of a non-overlapping concordant pair, 0 otherwise
	 */
	public final int pairMateStart;
	public CoverageRead(int referenceIndex, int alignmentStart, int alignmentEnd, int pairMateStart) {
		this.referenceIndex = referenceIndex;
		this.alignmentStart = alignmentStart;
		this.alignmentEnd = alignmentEnd;
		this.pairMateStart = pairMateStart;
	}
	/**
	 * Reads that count towards reference coverage: aligned and, if duplicates are ignored, not duplicates.
	 */
	public static boolean isCounted(SAMRecord read, boolean ignoreDuplicates) {
		return !read.getReadUnmappedFlag() && !(ignoreDuplicates && read.getDuplicateReadFlag());
	}
	public static CoverageRead create(SAMRecord read, ReadPairConcordanceCalculator pairing) {
		int pairMateStart = isLowerMappedOfNonOverlappingConcordantPair(read, pairing) ? read.getMateAlignmentStart() : 0;
		return new CoverageRead(read.getReferenceIndex(), read.getAlignmentStart(), read.getAlignmentEnd(), pairMateStart);
	}
	private static boolean isLowerMappedOfNonOverlappingConcordantPair(SAMRecord read, ReadPairConcordanceCalculator pairing) {
		return !read.getReadUnmappedFlag()
				&& read.getReadPairedFlag()
				&& !read.getMateUnmappedFlag()
				&& read.getAlignmentEnd() < read.getMateAlignmentStart()
				&& read.getReferenceIndex().equals(read.getMateReferenceIndex())
				&& (read.getAlignmentStart() < read.getMateAlignmentStart()
						|| (read.getAlignmentStart() == read.getMateAlignmentStart() && read.getFirstOfPairFlag()))
				&& pairing.isConcordant(read);
	}
}
