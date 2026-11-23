package gridss.analysis;

import htsjdk.samtools.CRAMFileWriter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.cram.structure.CRAMEncodingStrategy;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import org.junit.Assert;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class CollectTagMetricsSkippedTagsTest {
	private static final File REFERENCE = new File("src/test/resources/small.fa");
	@Rule
	public TemporaryFolder testFolder = new TemporaryFolder();

	private File writeCram() throws IOException {
		final SAMFileHeader header = new SAMFileHeader(SAMSequenceDictionaryExtractor.extractDictionary(REFERENCE.toPath()));
		header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
		header.addReadGroup(new SAMReadGroupRecord("rg"));
		final File cram = testFolder.newFile("tags.cram");
		try (FileOutputStream out = new FileOutputStream(cram)) {
			final CRAMFileWriter writer = new CRAMFileWriter(new CRAMEncodingStrategy(), out, null, true, new ReferenceSource(REFERENCE), header, cram.getName());
			for (int i = 0; i < 40; i++) {
				final SAMRecord r = new SAMRecord(header);
				r.setReadName("r" + i);
				r.setReferenceIndex(0);
				r.setAlignmentStart(1 + 10 * i);
				r.setCigarString("10M");
				r.setReadString("ACGTACGTAC");
				r.setBaseQualityString("IIIIIIIIII");
				r.setDuplicateReadFlag(i % 5 == 0);
				r.setAttribute("RG", "rg");
				if (i % 3 != 0) r.setAttribute("t0", "IIIII" + i);
				r.setAttribute("tp", new byte[] { (byte) i, 0, 1 });
				if (i % 2 == 0) r.setAttribute("NM", i);
				writer.addAlignment(r);
			}
			writer.close();
		}
		return cram;
	}

	private int runJava(final String skipTags, final File log, final String... args) throws IOException, InterruptedException {
		final List<String> cmd = new ArrayList<>(Arrays.asList(
				new File(System.getProperty("java.home"), "bin/java").getPath(),
				"-Dsamjdk.cram_skip_tags=" + skipTags, "-cp", System.getProperty("java.class.path")));
		cmd.addAll(Arrays.asList(args));
		return new ProcessBuilder(cmd).redirectErrorStream(true).redirectOutput(log).start().waitFor();
	}

	private File tagMetrics(final File cram, final String skipTags) throws IOException, InterruptedException {
		final File out = new File(testFolder.getRoot(), "skip_" + skipTags.replace(',', '_') + ".tag_metrics");
		Assert.assertEquals(0, runJava(skipTags, new File(out.getPath() + ".log"), CollectTagMetrics.class.getName(),
				"I=" + cram, "O=" + out, "R=" + REFERENCE.getAbsolutePath()));
		return out;
	}

	private static List<String> body(final File metrics) throws IOException {
		return Files.readAllLines(metrics.toPath()).stream().filter(l -> !l.startsWith("#")).collect(Collectors.toList());
	}

	@Test
	public void counts_skipped_tags_as_if_decoded() throws Exception {
		final File cram = writeCram();
		final File sam = new File(testFolder.getRoot(), "skipped.sam");
		Assert.assertEquals(0, runJava("t0,tp", sam, "picard.cmdline.PicardCommandLine", "ViewSam",
				"I=" + cram, "R=" + REFERENCE.getAbsolutePath(), "ALIGNMENT_STATUS=All", "PF_STATUS=All"));
		final String records = Files.readAllLines(sam.toPath()).stream().filter(l -> l.startsWith("r")).collect(Collectors.joining("\n"));
		Assert.assertFalse("t0 and tp must really be skipped", records.contains("t0:Z") || records.contains("tp:B"));
		Assert.assertTrue(records.contains("NM:i"));

		final List<String> expected = body(tagMetrics(cram, ""));
		Assert.assertTrue(expected.contains("t0\t21"));
		Assert.assertEquals(expected, body(tagMetrics(cram, "t0,tp")));
	}
}
