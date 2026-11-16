package gridss;

import au.edu.wehi.idsv.ProcessingContext;
import au.edu.wehi.idsv.SAMEvidenceSource;
import au.edu.wehi.idsv.SequentialCoverageAnnotator;
import au.edu.wehi.idsv.VariantContextDirectedEvidence;
import gridss.cmdline.VcfTransformCommandLineProgram;
import htsjdk.samtools.util.CloseableIterator;

import java.util.List;
import java.util.concurrent.ExecutorService;

public class AnnotateReferenceCoverage extends VcfTransformCommandLineProgram {
	/**
	 * Defensive programming safety margin around expected window size
	 */
	private static final int WINDOW_SIZE_SAFETY_MARGIN = 100000;
	public static int coverageWindowSize(ProcessingContext context, List<SAMEvidenceSource> sources) {
		return 2 * SAMEvidenceSource.maximumWindowSize(context, sources, null) + WINDOW_SIZE_SAFETY_MARGIN;
	}
	@Override
	public CloseableIterator<VariantContextDirectedEvidence> iterator(CloseableIterator<VariantContextDirectedEvidence> calls, ExecutorService threadpool) {
		ProcessingContext context = getContext();
		List<SAMEvidenceSource> sources = getSamEvidenceSources();
		return new SequentialCoverageAnnotator<VariantContextDirectedEvidence>(context, sources, calls, coverageWindowSize(context, sources), threadpool);
	}
	public static void main(String[] argv) {
        System.exit(new AnnotateReferenceCoverage().instanceMain(argv));
    }
}
