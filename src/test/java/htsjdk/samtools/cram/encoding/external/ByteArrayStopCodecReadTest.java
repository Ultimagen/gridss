package htsjdk.samtools.cram.encoding.external;

import org.junit.Assert;
import org.junit.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

public class ByteArrayStopCodecReadTest {
    private static final byte STOP = '\t';

    @Test
    public void readsBackWrittenArraysOfAnyLength() {
        final Random random = new Random(3);
        final List<byte[]> values = new ArrayList<>();
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        final ByteArrayStopCodec writer = new ByteArrayStopCodec(null, out, STOP);
        for (int i = 0; i < 500; i++) {
            final byte[] v = new byte[random.nextInt(i % 50 == 0 ? 3000 : 300)];
            for (int j = 0; j < v.length; j++) {
                v[j] = (byte) ('A' + random.nextInt(26));
            }
            values.add(v);
            writer.write(v);
        }
        final byte[] trailing = {1, 2, 3};
        out.write(trailing, 0, trailing.length);
        final ByteArrayInputStream in = new ByteArrayInputStream(out.toByteArray());
        final ByteArrayStopCodec reader = new ByteArrayStopCodec(in, null, STOP);
        for (final byte[] v : values) {
            Assert.assertArrayEquals(v, reader.read());
        }
        Assert.assertArrayEquals(trailing, reader.read());
        Assert.assertArrayEquals(new byte[0], reader.read());
        Assert.assertEquals(-1, in.read());
    }

    @Test
    public void leavesStreamJustAfterTheStopByte() {
        final ByteArrayInputStream in = new ByteArrayInputStream(new byte[]{'a', 'b', STOP, 'c', STOP, STOP, 'd'});
        final ByteArrayStopCodec reader = new ByteArrayStopCodec(in, null, STOP);
        Assert.assertArrayEquals("ab".getBytes(), reader.read());
        Assert.assertEquals('c', in.read());
        Assert.assertEquals(STOP, in.read());
        Assert.assertArrayEquals(new byte[0], reader.read());
        Assert.assertArrayEquals(new byte[]{'d'}, reader.read());
        Assert.assertTrue(Arrays.equals(new byte[0], reader.read()));
    }
}
