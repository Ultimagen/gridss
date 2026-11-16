package htsjdk.samtools.cram.compression.rans;

import org.junit.Assert;
import org.junit.Test;

import java.nio.ByteBuffer;
import java.util.Random;

public class RANSDecoderResetTest {
    private static byte[] data(final Random random, final int length, final int alphabet) {
        final byte[] b = new byte[length];
        for (int i = 0; i < length; i++) {
            b[i] = (byte) (i > 0 && random.nextInt(4) == 0 ? b[i - 1] : 'A' + random.nextInt(alphabet));
        }
        return b;
    }

    private static byte[] bytes(final ByteBuffer buffer) {
        final byte[] b = new byte[buffer.remaining()];
        buffer.get(b);
        return b;
    }

    @Test
    public void sharedInstanceDecodesInterleavedOrdersAndAlphabets() {
        final Random random = new Random(42);
        final RANS encoder = new RANS();
        final RANS decoder = new RANS();
        for (int round = 0; round < 200; round++) {
            final byte[] in = data(random, 1 + random.nextInt(5000), 1 + random.nextInt(round % 3 == 0 ? 4 : 60));
            final RANS.ORDER order = random.nextBoolean() ? RANS.ORDER.ZERO : RANS.ORDER.ONE;
            final ByteBuffer compressed = encoder.compress(ByteBuffer.wrap(in), order);
            Assert.assertArrayEquals("round " + round, in, bytes(decoder.uncompress(compressed)));
        }
    }

    @Test
    public void encodeAfterDecodeOnOneInstance() {
        final Random random = new Random(7);
        final RANS rans = new RANS();
        for (int round = 0; round < 50; round++) {
            final byte[] in = data(random, 100 + random.nextInt(2000), 20);
            final RANS.ORDER order = round % 2 == 0 ? RANS.ORDER.ZERO : RANS.ORDER.ONE;
            final ByteBuffer compressed = rans.compress(ByteBuffer.wrap(in), order);
            Assert.assertArrayEquals(in, bytes(rans.uncompress(compressed)));
        }
    }
}
