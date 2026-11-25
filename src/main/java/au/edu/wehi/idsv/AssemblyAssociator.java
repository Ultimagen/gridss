package au.edu.wehi.idsv;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.PriorityQueue;

/**
 * Maps reads to their associated assemblies.
 * 
 * Note: this iterator requires that the exact evidence supplied to the
 * original assemblies is supplied to this iterator. Evidence IDs that never
 * arrive are dropped once their assembly is more than windowSize behind.
 * 
 * @author Daniel Cameron
 *
 */
public class AssemblyAssociator implements CloseableIterator<DirectedEvidence> {
	private final Iterator<DirectedEvidence> it;
	private final Iterator<SAMRecord> assit;
	private final int windowSize;
	private final HashMap<String, String> evidenceToAssemblyName = new HashMap<>();
	private final PriorityQueue<LoadedAssembly> loaded = new PriorityQueue<>(
			Comparator.<LoadedAssembly>comparingInt(a -> a.referenceIndex).thenComparingInt(a -> a.unclippedEnd));
	private SAMRecord lastAssembly = null;
	public AssemblyAssociator(Iterator<DirectedEvidence> it, Iterator<SAMRecord> rawAssemblies, int windowSize) {
		this.it = it;
		this.assit = rawAssemblies;
		this.windowSize = windowSize;
	}
	@Override
	public boolean hasNext() {
		return it.hasNext();
	}

	@Override
	public DirectedEvidence next() {
		DirectedEvidence e = it.next();
		associate(e);
		return e;
	}
	private DirectedEvidence associate(DirectedEvidence e) {
		if (e == null || e.getBreakendSummary() == null) {
			return e;
		}
		ensureAssembliesLoadedUntil(e.getBreakendSummary());
		setAssociatedAssembly(e, evidenceToAssemblyName.remove(e.getEvidenceID()));
		flushBefore(e.getBreakendSummary());
		return e;
	}
	/**
	 * Drops the evidence IDs of assemblies more than windowSize behind, the mirror of the loading window.
	 * Externally built assemblies (sv_rematch) list reads that never arrive here, so their IDs would otherwise stay until close.
	 */
	private void flushBefore(BreakendSummary breakendSummary) {
		while (!loaded.isEmpty() && isBefore(loaded.peek(), breakendSummary)) {
			LoadedAssembly a = loaded.poll();
			for (String eid : a.evidenceIDs) {
				evidenceToAssemblyName.remove(eid, a.name);
			}
		}
	}
	private boolean isBefore(LoadedAssembly a, BreakendSummary breakendSummary) {
		return a.referenceIndex < breakendSummary.referenceIndex ||
			(a.referenceIndex == breakendSummary.referenceIndex && a.unclippedEnd + windowSize < breakendSummary.start);
	}
	private void ensureAssembliesLoadedUntil(BreakendSummary breakendSummary) {
		assert(breakendSummary != null);
		while (assit.hasNext() && (lastAssembly == null || !isAfter(breakendSummary, lastAssembly))) {
			lastAssembly = assit.next();
			load(lastAssembly);
		}
	}
	private void load(SAMRecord ass) {
		assert(ass != null);
		Collection<String> eids = new AssemblyAttributes(ass).getEvidenceIDs(null, null, null, null);
		for (String eid : eids) {
			evidenceToAssemblyName.put(eid, ass.getReadName());
		}
		if (!ass.getReadUnmappedFlag()) {
			loaded.add(new LoadedAssembly(ass.getReferenceIndex(), ass.getUnclippedEnd(), ass.getReadName(), eids));
		}
	}
	private static class LoadedAssembly {
		final int referenceIndex;
		final int unclippedEnd;
		final String name;
		final Collection<String> evidenceIDs;
		LoadedAssembly(int referenceIndex, int unclippedEnd, String name, Collection<String> evidenceIDs) {
			this.referenceIndex = referenceIndex;
			this.unclippedEnd = unclippedEnd;
			this.name = name;
			this.evidenceIDs = evidenceIDs;
		}
	}
	private boolean isAfter(BreakendSummary breakendSummary, SAMRecord position) {
		return position.getReferenceIndex() > breakendSummary.referenceIndex ||
			(position.getReferenceIndex() == breakendSummary.referenceIndex && position.getUnclippedStart() > breakendSummary.end + windowSize);
	}
	@Override
	public void close() {
		CloserUtil.close(it);
		CloserUtil.close(assit);
	}
	private void setAssociatedAssembly(DirectedEvidence e, String assemblyName) {
		if (e instanceof SingleReadEvidence) {
			((SingleReadEvidence)e).setAssociatedAssemblyName(assemblyName);
		} else if (e instanceof NonReferenceReadPair) {
			((NonReferenceReadPair)e).setAssociatedAssemblyName(assemblyName);
		}
	}
}
