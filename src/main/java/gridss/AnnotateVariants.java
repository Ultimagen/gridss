package gridss;

import au.edu.wehi.idsv.*;
import com.google.common.io.Files;
import gridss.cmdline.VcfTransformCommandLineProgram;
import htsjdk.samtools.util.AsyncBufferedIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "Annotates breakpoint variant calls",  
        oneLineSummary = "Annotates breakpoint variant calls",
        programGroup = gridss.cmdline.programgroups.VariantCalling.class
)
public class AnnotateVariants extends VcfTransformCommandLineProgram {
	private static final Log log = Log.getInstance(AnnotateVariants.class);
	/**
	 * Decode each input once for both evidence allocation and reference coverage
	 */
	public boolean sharedCoverageDecode = Defaults.ANNOTATE_SHARED_DECODE;
	public int sharedDecodeMaxCoverageReads = Defaults.ANNOTATE_SHARED_DECODE_MAX_COVERAGE_READS;
	public int sharedDecodeMaxEvidenceRecords = Defaults.ANNOTATE_SHARED_DECODE_MAX_EVIDENCE_RECORDS;
	public static void writeAssemblyBreakends(File file, AssemblyEvidenceSource assemblyEvidence) throws IOException {
		log.info("Writing breakend assembly support.");
		File tmp = gridss.Defaults.OUTPUT_TO_TEMP_FILE ? FileSystemContext.getWorkingFileFor(file) : file;
		try (BufferedOutputStream writer = new BufferedOutputStream(new FileOutputStream(tmp))) {
			try (CloseableIterator<DirectedEvidence> it = assemblyEvidence.iterator(SAMEvidenceSource.EvidenceSortOrder.EvidenceStartPosition)) {
				while (it.hasNext()) {
					SingleReadEvidence ass = (SingleReadEvidence) it.next();
					if (!ass.getSAMRecord().isSecondaryOrSupplementary()) {
						writer.write('>');
						writer.write(ass.getEvidenceID().getBytes(StandardCharsets.US_ASCII));
						writer.write('\n');
						if (ass.getBreakendSummary().direction == BreakendDirection.Forward) {
							writer.write(ass.getAnchorSequence(), 0, ass.getAnchorSequence().length);
							writer.write(ass.getBreakendSequence(), 0, ass.getBreakendSequence().length);
						} else {
							writer.write(ass.getBreakendSequence(), 0, ass.getBreakendSequence().length);
							writer.write(ass.getAnchorSequence(), 0, ass.getAnchorSequence().length);
						}
						writer.write('\n');
					}
				}
			}
		}
		if (tmp != file) {
			Files.move(tmp, file);
		}
		log.info("Writing breakend assembly support complete.");
	}
	@Override
	public CloseableIterator<VariantContextDirectedEvidence> iterator(CloseableIterator<VariantContextDirectedEvidence> calls, ExecutorService threadpool) {
		AllocateEvidence ae = new AllocateEvidence();
		AnnotateReferenceCoverage arc = new AnnotateReferenceCoverage();
		AnnotateInexactHomology ihom = new AnnotateInexactHomology();
		copyInputs(ae);
		copyInputs(arc);
		copyInputs(ihom);
		ae.INPUT_VCF = INPUT_VCF; // needed for caching 
		if (sharedCoverageDecode) {
			shareDecodeWithCoverage();
		}
		calls = new AsyncBufferedIterator<VariantContextDirectedEvidence>(ae.iterator(calls, threadpool), 128, 2, "AllocateEvidence");
		//List<VariantContextDirectedEvidence> callsList = calls.stream().collect(Collectors.toList());
		calls = new AsyncBufferedIterator<VariantContextDirectedEvidence>(arc.iterator(calls, threadpool), 128, 2, "AnnotateReferenceCoverage");

		calls = ihom.iterator(calls, threadpool);
		AnnotateAlleleFraction aaf = new AnnotateAlleleFraction(new AlleleFractionAnnotator(getContext(), getSamEvidenceSources()));
		return aaf.iterator(calls, threadpool);
	}
	/**
	 * Must run before ae.iterator(): its read-ahead threads start decoding straight away.
	 */
	private void shareDecodeWithCoverage() {
		int windowSize = AnnotateReferenceCoverage.coverageWindowSize(getContext(), getSamEvidenceSources());
		// allocation can move a call within its evidence window; windowSize is far larger
		List<int[]> queries = coverageQueries(windowSize);
		for (SAMEvidenceSource ses : getSamEvidenceSources()) {
			CoverageReadTee tee = ses.shareDecodeWithCoverage(sharedDecodeMaxCoverageReads, sharedDecodeMaxEvidenceRecords);
			if (tee == null) {
				log.info("Not sharing the decode of " + ses.getFile() + " with reference coverage: evidence reads a different file");
				continue;
			}
			int maxEvidenceWindow = SequentialReferenceCoverageLookup.maxEvidenceWindow(ses.getMetrics().getIdsvMetrics(), ses.getReadPairConcordanceCalculator());
			tee.setCountableWindows(new CoverageReadWindows(queries, windowSize, maxEvidenceWindow, getContext().getDictionary().size()));
		}
	}
	private List<int[]> coverageQueries(int margin) {
		List<int[]> queries = new ArrayList<>();
		try (CloseableIterator<VariantContextDirectedEvidence> it = getBreakends(INPUT_VCF)) {
			while (it.hasNext()) {
				BreakendSummary bs = it.next().getBreakendSummary();
				queries.add(new int[] { bs.referenceIndex, bs.start - 1 - margin, bs.end + margin });
			}
		}
		return queries;
	}
	public static void main(String[] argv) {
        System.exit(new AnnotateVariants().instanceMain(argv));
    }
}
