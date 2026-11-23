package htsjdk.samtools.cram.encoding.external;

import htsjdk.samtools.util.RuntimeIOException;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.Arrays;

/**
 * Encode byte arrays by specifying a stop byte to separate the arrays.
 * This cannot be a byte that appears in the data.
 */
final class ByteArrayStopCodec extends ExternalCodec<byte[]> {
    private final int stop;
    private byte[] scratch = new byte[256];

    /**
     * Construct a Byte Array Stop Codec
     *
     * @param inputStream the input bytestream to read from
     * @param outputStream the output bytestream to write to
     * @param stopByte the byte used to mark array boundaries
     */
    public ByteArrayStopCodec(final ByteArrayInputStream inputStream,
                              final ByteArrayOutputStream outputStream,
                              final byte stopByte) {
        super(inputStream, outputStream);
        this.stop = 0xFF & stopByte;
    }

    @Override
    public byte[] read() {
        // Ultima fork: bulk read + scan instead of one synchronized read() per byte (~25% of CRAM decode CPU).
        // Same result and stream position: the bytes before the stop byte, with the stop byte consumed.
        inputStream.mark(0);
        int n = 0;
        while (true) {
            if (n == scratch.length) {
                scratch = Arrays.copyOf(scratch, 2 * n);
            }
            final int r = inputStream.read(scratch, n, scratch.length - n);
            if (r <= 0) {
                return Arrays.copyOf(scratch, n);
            }
            for (int i = n; i < n + r; i++) {
                if ((0xFF & scratch[i]) == stop) {
                    inputStream.reset();
                    inputStream.skip(i + 1);
                    return Arrays.copyOf(scratch, i);
                }
            }
            n += r;
        }
    }

    @Override
    public byte[] read(final int length) {
        throw new RuntimeException("Not implemented.");
    }

    @Override
    public void write(final byte[] value) {
        try {
            outputStream.write(value);
            outputStream.write(stop);
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }
}
