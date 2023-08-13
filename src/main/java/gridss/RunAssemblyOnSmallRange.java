package gridss;
import au.edu.wehi.idsv.*;
import au.edu.wehi.idsv.bed.IntervalBed;
import au.edu.wehi.idsv.debruijn.positional.PositionalAssembler;
import gridss.cmdline.MultipleSamFileCommandLineProgram;
import htsjdk.samtools.BAMFileWriter;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.util.concurrent.ExecutorService;

public class RunAssemblyOnSmallRange extends MultipleSamFileCommandLineProgram {
    private static final Log log = Log.getInstance(CallVariants.class);

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output file containing subset of input", optional = false)
    public File OUTPUT;
    public static final File CONFIGURATION_FILE = new File("gridss.config");
    public static void main(String[] argv) {
        System.exit(new RunAssemblyOnSmallRange().instanceMain(argv));
    }

    public int doWork() {
        IOUtil.assertFileIsWritable(OUTPUT);
        ProcessingContext pc = getContext();
        IntervalBed excluded = new IntervalBed(pc.getLinear());
        AssemblyEvidenceSource aes = new AssemblyEvidenceSource(pc, getSamEvidenceSources(), OUTPUT);
        AggregateEvidenceSource des = new AggregateEvidenceSource(pc, getSamEvidenceSources(),null, SAMEvidenceSource.EvidenceSortOrder.EvidenceStartPosition);

        PositionalAssembler pose = new PositionalAssembler(pc, aes, new SequentialIdGenerator("asm"), des.iterator(), BreakendDirection.Forward, excluded, null);
        File outputFile = new File(OUTPUT,"w");
        SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(aes.getHeader(), false)) {;
        while (pose.hasNext()) {
            SAMRecord asm = pose.next();
            log.info(asm);
        }
        return 0;
    }

    public int doWork(ExecutorService threadpool){
        return 0;
    }
}
