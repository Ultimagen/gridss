package au.edu.wehi.idsv;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Read start positions that a {@link SequentialReferenceCoverageLookup} can count for a known set of queries.
 *
 * The lookup skips, without counting, every read starting more than its window plus its maximum
 * evidence window before the position it is asked about. A read that no query can reach that way
 * is dropped before it is queued, so the queue between a shared decode and the coverage lookup
 * holds only reads near the calls being annotated.
 *
 * Lookups must be in coordinate order.
 */
public class CoverageReadWindows {
	private final int[][] starts;
	private final int[][] ends;
	private int currentReferenceIndex = -1;
	private int index = 0;
	/**
	 * @param queries per call, the reference index and the first and last position the lookup may be asked about
	 * @param windowSize lookup window size
	 * @param maxEvidenceWindow lookup maximum evidence window
	 * @param referenceCount number of reference sequences
	 */
	public CoverageReadWindows(List<int[]> queries, int windowSize, int maxEvidenceWindow, int referenceCount) {
		List<List<int[]>> byReference = new ArrayList<>();
		for (int i = 0; i < referenceCount; i++) byReference.add(new ArrayList<>());
		for (int[] q : queries) {
			long start = (long) q[1] - windowSize - maxEvidenceWindow - 1;
			byReference.get(q[0]).add(new int[] { (int) Math.max(Integer.MIN_VALUE, start), q[2] });
		}
		starts = new int[referenceCount][];
		ends = new int[referenceCount][];
		for (int i = 0; i < referenceCount; i++) {
			List<int[]> list = byReference.get(i);
			list.sort((a, b) -> Integer.compare(a[0], b[0]));
			int[] s = new int[list.size()];
			int[] e = new int[list.size()];
			int n = 0;
			for (int[] w : list) {
				if (n > 0 && w[0] <= e[n - 1] + 1) {
					e[n - 1] = Math.max(e[n - 1], w[1]);
				} else {
					s[n] = w[0];
					e[n] = w[1];
					n++;
				}
			}
			starts[i] = Arrays.copyOf(s, n);
			ends[i] = Arrays.copyOf(e, n);
		}
	}
	public boolean isCountable(int referenceIndex, int alignmentStart) {
		if (referenceIndex < 0 || referenceIndex >= starts.length) return false;
		if (referenceIndex != currentReferenceIndex) {
			currentReferenceIndex = referenceIndex;
			index = 0;
		}
		int[] s = starts[referenceIndex];
		int[] e = ends[referenceIndex];
		while (index < e.length && e[index] < alignmentStart) index++;
		return index < s.length && s[index] <= alignmentStart;
	}
}
