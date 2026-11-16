package gridss;

import au.edu.wehi.idsv.*;
import au.edu.wehi.idsv.util.AutoClosingIterator;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.util.concurrent.MoreExecutors;
import htsjdk.samtools.SAMRecord;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.stream.Collectors;

import static org.junit.Assert.*;

public class AnnotateVariantsSharedDecodeTest extends IntermediateFilesTest {
	private File assemblyFile;
	private ProcessingContext pc;
	private void setupInput() throws IOException {
		pc = getCommandlineContext();
		pc.getVariantCallingParameters().minSize = 0;
		pc.getVariantCallingParameters().minScore = 0;
		pc.getVariantCallingParameters().minReads = 0;
		pc.getVariantCallingParameters().writeFiltered = true;
		Random rng = new Random(0);
		List<SAMRecord> reads = new ArrayList<>();
		for (int i = 0; i < 3000; i++) {
			int pos = 1 + rng.nextInt(2000);
			for (SAMRecord r : RP(0, pos, pos + 20 + rng.nextInt(40), 10)) {
				r.setReadName(r.getReadName() + "#" + i);
				r.setDuplicateReadFlag(rng.nextInt(20) == 0);
				reads.add(r);
			}
		}
		for (int i = 0; i < 40; i++) {
			int pos = 500 + rng.nextInt(3);
			for (SAMRecord r : DP(0, pos, "5M5S", true, 1, 100 + rng.nextInt(5), "5M", true)) {
				r.setReadName(r.getReadName() + "#dp" + i);
				reads.add(r);
			}
			for (SAMRecord r : DP(0, 1500 + rng.nextInt(3), "10M", true, 0, 1800 + rng.nextInt(3), "10M", true)) {
				r.setReadName(r.getReadName() + "#dp2" + i);
				reads.add(r);
			}
		}
		createInput(reads);
		SAMEvidenceSource ses = new SAMEvidenceSource(getContext(), input, null, 0);
		ses.ensureMetrics();
		assemblyFile = new File(testFolder.getRoot(), "assembly.bam");
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), assemblyFile);
		aes.assembleBreakends(null);
		aes.ensureExtracted();
		VariantCaller caller = new VariantCaller(pc, ImmutableList.of(ses), ImmutableList.of(aes));
		caller.callBreakends(output, MoreExecutors.newDirectExecutorService());
	}
	private List<String> annotate(boolean shared, int maxCoverageReads, int maxEvidenceRecords, List<SAMEvidenceSource> sources) {
		SAMEvidenceSource ses = new SAMEvidenceSource(getContext(), input, null, 0);
		ses.ensureMetrics();
		sources.add(ses);
		AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, ImmutableList.of(ses), assemblyFile);
		AnnotateVariants cmd = new AnnotateVariants();
		cmd.sharedCoverageDecode = shared;
		cmd.sharedDecodeMaxCoverageReads = maxCoverageReads;
		cmd.sharedDecodeMaxEvidenceRecords = maxEvidenceRecords;
		cmd.INPUT_VCF = output;
		cmd.setContext(pc);
		cmd.setAssemblySource(ImmutableList.of(aes));
		cmd.setSamEvidenceSources(ImmutableList.of(ses));
		List<VariantContextDirectedEvidence> calls = Lists.newArrayList(Iterables.filter(getVcf(output, null), VariantContextDirectedEvidence.class));
		ExecutorService threadpool = Executors.newFixedThreadPool(2);
		try {
			return Lists.newArrayList(cmd.iterator(new AutoClosingIterator<>(calls.iterator()), threadpool))
					.stream()
					.map(v -> v.toStringDecodeGenotypes().replaceAll("\\[I@[0-9a-f]+", "")
							+ v.getAttributeAsIntList("CIPOS", 0) + v.getAttributeAsIntList("CIRPOS", 0))
					.collect(Collectors.toList());
		} finally {
			threadpool.shutdown();
		}
	}
	@Test
	public void shared_decode_should_give_the_same_annotations() throws IOException {
		setupInput();
		List<SAMEvidenceSource> sources = new ArrayList<>();
		List<String> expected = annotate(false, 0, 0, sources);
		assertTrue(expected.size() > 4);
		assertTrue(expected.stream().anyMatch(s -> s.matches(".*attr=\\{.* REF=[1-9].*")));
		assertTrue(expected.stream().anyMatch(s -> s.matches(".*attr=\\{.* REFPAIR=[1-9].*")));
		assertNull(sources.get(0).getCoverageReadTee());

		assertEquals(expected, annotate(true, 1 << 20, 10000, sources));
		CoverageReadTee tee = sources.get(1).getCoverageReadTee();
		assertNull("evidence used the shared decode", tee.claimEvidenceIterator());
		assertNull("coverage used the shared decode", tee.claimCoverageIterator());
		assertEquals(-1, tee.getDetachedAt());

		assertEquals(expected, annotate(true, 50, 10000, sources));
		assertTrue(sources.get(2).getCoverageReadTee().getDetachedAt() > 0);
	}
}
