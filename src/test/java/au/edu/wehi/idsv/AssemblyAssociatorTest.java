package au.edu.wehi.idsv;

import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMRecord;
import org.junit.Test;

import java.lang.reflect.Field;
import java.util.Map;
import java.util.List;

import static org.junit.Assert.assertEquals;

public class AssemblyAssociatorTest extends TestHelper {
	private static final int WINDOW = 10;
	private final SequentialIdGenerator ids = new SequentialIdGenerator("asm");

	private SoftClipEvidence evidence(int referenceIndex, int position) {
		return SCE(FWD, Read(referenceIndex, position, "1M2S"));
	}
	private SoftClipEvidence otherRead(int referenceIndex, int position) {
		return SCE(FWD, Read(referenceIndex, position, "1M3S"));
	}
	private SAMRecord assembly(int referenceIndex, int position, DirectedEvidence... support) {
		return AssemblyFactory.createUnanchoredBreakend(getContext(), AES(), ids,
				new BreakendSummary(referenceIndex, FWD, position, position, position),
				ImmutableList.copyOf(support), null, B("GTAC"), new byte[] {1, 2, 3, 4});
	}
	private AssemblyAssociator associator(List<DirectedEvidence> evidence, SAMRecord... assemblies) {
		return new AssemblyAssociator(evidence.iterator(), ImmutableList.copyOf(assemblies).iterator(), WINDOW);
	}

	private static int retained(AssemblyAssociator it) throws ReflectiveOperationException {
		Field f = AssemblyAssociator.class.getDeclaredField("evidenceToAssemblyName");
		f.setAccessible(true);
		return ((Map<?, ?>) f.get(it)).size();
	}

	@Test
	public void should_associate_evidence_with_its_assembly() {
		SoftClipEvidence e = evidence(0, 5);
		SAMRecord a = assembly(0, 5, e);
		AssemblyAssociator it = associator(ImmutableList.of(e), a);
		it.next();
		assertEquals(a.getReadName(), e.getAssociatedAssemblyName());
	}

	@Test
	public void should_drop_evidence_ids_that_never_arrive_once_the_assembly_is_window_behind() throws ReflectiveOperationException {
		SoftClipEvidence neverSupplied = otherRead(0, 5);
		SoftClipEvidence near = evidence(0, 5);
		SoftClipEvidence far = evidence(0, 5 + 3 * WINDOW);
		AssemblyAssociator it = associator(ImmutableList.of(near, far), assembly(0, 5, neverSupplied, near));
		it.next();
		assertEquals(1, retained(it));
		it.next();
		assertEquals(0, retained(it));
	}

	@Test
	public void should_keep_evidence_ids_while_the_assembly_is_within_the_window() throws ReflectiveOperationException {
		SoftClipEvidence neverSupplied = otherRead(0, 5);
		SoftClipEvidence near = evidence(0, 5);
		SoftClipEvidence stillNear = evidence(0, 5 + WINDOW / 2);
		AssemblyAssociator it = associator(ImmutableList.of(near, stillNear), assembly(0, 5, neverSupplied, near));
		it.next();
		it.next();
		assertEquals(1, retained(it));
	}

	@Test
	public void should_drop_evidence_ids_of_an_earlier_contig() throws ReflectiveOperationException {
		SoftClipEvidence neverSupplied = otherRead(0, 5);
		SoftClipEvidence nextContig = evidence(1, 1);
		AssemblyAssociator it = associator(ImmutableList.of(nextContig), assembly(0, 5, neverSupplied));
		it.next();
		assertEquals(0, retained(it));
	}

	@Test
	public void dropping_an_older_assembly_should_keep_a_shared_read_on_the_newer_one() {
		SoftClipEvidence shared = evidence(0, 4 + 3 * WINDOW);
		SAMRecord older = assembly(0, 5, shared);
		SAMRecord newer = assembly(0, 4 + 3 * WINDOW, shared);
		SoftClipEvidence beforeShared = evidence(0, 3 + 3 * WINDOW);
		AssemblyAssociator it = associator(ImmutableList.of(beforeShared, shared), older, newer);
		it.next();
		it.next();
		assertEquals(newer.getReadName(), shared.getAssociatedAssemblyName());
	}
}
