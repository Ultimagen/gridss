package htsjdk.samtools.cram.encoding.reader;

import htsjdk.samtools.cram.io.ITF8;
import htsjdk.samtools.cram.structure.CompressionHeader;
import htsjdk.samtools.cram.structure.DataSeries;
import htsjdk.samtools.cram.structure.EncodingDescriptor;
import htsjdk.samtools.cram.structure.EncodingID;
import htsjdk.samtools.cram.structure.ReadTag;

import java.nio.ByteBuffer;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * GRIDSS fork: the tags of a container that CRAM decode can leave out (samjdk.cram_skip_tags).
 *
 * A tag is skipped only if all its data sits in external blocks that no other data series or kept tag reads,
 * so leaving those blocks compressed and unread changes nothing else in the decode.
 */
public final class SkippedTags {
    /**
     * Transient SAMRecord attribute: the String[] names of the tags decode left out of that record, for readers
     * that only count tags (CollectTagMetrics). Shared between records, don't modify.
     */
    public static final String SKIPPED_TAGS_ATTRIBUTE = "samjdk.cram_skipped_tags";

    /** Tag IDs as in {@link htsjdk.samtools.cram.structure.ReadTag#name3BytesToInt(byte[])} */
    final Set<Integer> tagIds = new HashSet<>();
    final Set<Integer> contentIds = new HashSet<>();

    SkippedTags(final CompressionHeader header, final Set<String> names) {
        if (names.isEmpty()) {
            return;
        }
        final Map<Integer, Set<Integer>> candidates = new HashMap<>();
        final Set<Integer> readByOthers = new HashSet<>();
        for (final Map.Entry<Integer, EncodingDescriptor> e : header.getTagEncodingMap().entrySet()) {
            final Set<Integer> ids = new HashSet<>();
            final boolean externalOnly = collectExternalIds(e.getValue().getEncodingID(), e.getValue().getEncodingParameters(), ids);
            if (externalOnly && names.contains(tagName(e.getKey()))) {
                candidates.put(e.getKey(), ids);
            } else {
                readByOthers.addAll(ids);
            }
        }
        for (final DataSeries dataSeries : DataSeries.values()) {
            final EncodingDescriptor d = header.getEncodingMap().getEncodingDescriptorForDataSeries(dataSeries);
            if (d != null) {
                collectExternalIds(d.getEncodingID(), d.getEncodingParameters(), readByOthers);
            }
        }
        for (final Map.Entry<Integer, Set<Integer>> c : candidates.entrySet()) {
            if (c.getValue().stream().noneMatch(readByOthers::contains)) {
                tagIds.add(c.getKey());
                contentIds.addAll(c.getValue());
            }
        }
    }

    /**
     * @return for each tag ID list of the header's dictionary, the names of its skipped tags; null if none is skipped
     */
    public static String[][] namesByTagIdList(final CompressionHeader header, final Set<String> names) {
        final Set<Integer> skipped = new SkippedTags(header, names).tagIds;
        if (skipped.isEmpty()) {
            return null;
        }
        final byte[][][] dictionary = header.getTagIDDictionary();
        final String[][] byList = new String[dictionary.length][];
        for (int i = 0; i < dictionary.length; i++) {
            byList[i] = Arrays.stream(dictionary[i])
                    .map(ReadTag::name3BytesToInt)
                    .filter(skipped::contains)
                    .map(SkippedTags::tagName)
                    .toArray(String[]::new);
        }
        return byList;
    }

    static String tagName(final int tagId) {
        return new String(new char[] { (char) ((tagId >> 16) & 0xFF), (char) ((tagId >> 8) & 0xFF) });
    }

    /**
     * Adds the external blocks the encoding reads to ids
     * @return true if it reads no other block
     */
    private static boolean collectExternalIds(final EncodingID id, final byte[] params, final Set<Integer> ids) {
        final ByteBuffer buffer = ByteBuffer.wrap(params);
        switch (id) {
            case EXTERNAL:
                ids.add(ITF8.readUnsignedITF8(buffer));
                return true;
            case BYTE_ARRAY_STOP:
                buffer.get(); // stop byte
                ids.add(ITF8.readUnsignedITF8(buffer));
                return true;
            case BYTE_ARRAY_LEN:
                // as in ByteArrayLenEncoding.fromSerializedEncodingParams; both sub-encodings are always read
                final boolean lengths = collectSubEncoding(buffer, ids);
                final boolean values = collectSubEncoding(buffer, ids);
                return lengths && values;
            case HUFFMAN:
                return readsNoBits(buffer);
            default:
                return false;
        }
    }

    /**
     * A one-symbol Huffman code has bit length 0: it never reads the core block (e.g. htsjdk's fixed array lengths)
     */
    private static boolean readsNoBits(final ByteBuffer buffer) {
        final int symbols = ITF8.readUnsignedITF8(buffer);
        for (int i = 0; i < symbols; i++) {
            ITF8.readUnsignedITF8(buffer);
        }
        final int lengths = ITF8.readUnsignedITF8(buffer);
        for (int i = 0; i < lengths; i++) {
            if (ITF8.readUnsignedITF8(buffer) != 0) {
                return false;
            }
        }
        return true;
    }

    private static boolean collectSubEncoding(final ByteBuffer buffer, final Set<Integer> ids) {
        final EncodingID id = EncodingID.values()[buffer.get()];
        final byte[] params = new byte[ITF8.readUnsignedITF8(buffer)];
        buffer.get(params);
        return collectExternalIds(id, params, ids);
    }
}
