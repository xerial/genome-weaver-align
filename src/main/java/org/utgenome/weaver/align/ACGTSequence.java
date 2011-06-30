/*--------------------------------------------------------------------------
 *  Copyright 2011 utgenome.org
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *--------------------------------------------------------------------------*/
//--------------------------------------
// genome-weaver Project
//
// ACGTSequence.java
// Since: 2011/06/30
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

import org.utgenome.gwt.utgb.client.bio.IUPAC;
import org.xerial.snappy.SnappyInputStream;
import org.xerial.snappy.SnappyOutputStream;
import org.xerial.util.log.Logger;

/**
 * ACGT genome sequence with N, using 3-bit encoding
 * 
 * @author leo
 * 
 */
public class ACGTSequence implements LSeq
{
    private static Logger     _logger        = Logger.getLogger(ACGTSequence.class);

    private static final int  LONG_BYTE_SIZE = 8;
    private static final long MAX_SIZE       = (2 * 1024 * 1024 * 1024 * 8 * 8) / 3;

    // 2G (max of Java array size) * 8 (long byte size) * 8 / 3 (bit) =  42G  (42G characters)
    private long[]            seq;
    private long              numBases;

    public ACGTSequence(String s) {
        this(s.length());

        for (int i = 0; i < s.length(); ++i) {
            byte code = ACGT.to3bitCode(s.charAt(i));
            set(i, code);
        }
    }

    public ACGTSequence(long numBases) {
        this.numBases = numBases;
        long bitSize = numBases * 3;
        long blockBitSize = LONG_BYTE_SIZE * 3 * 8;
        long arraySize = ((bitSize + blockBitSize - 1) / blockBitSize) * 3;
        if (arraySize > Integer.MAX_VALUE) {
            throw new IllegalArgumentException(String.format("Cannot create ACGTSequece more than %,d size: %,d",
                    MAX_SIZE, numBases));
        }

        seq = new long[(int) arraySize];
    }

    private ACGTSequence(long[] rawSeq, long numBases) {
        this.seq = rawSeq;
        this.numBases = numBases;
    }

    public ACGT getACGT(long index) {
        return ACGT.decode((byte) lookup(index));
    }

    @Override
    public long lookup(long index) {
        int pos = (int) (index >> 6);
        int offset = (int) (index & 0x03FL);
        int shift = 62 - ((int) (index & 0x1FL) << 1);

        int nFlag = (int) (seq[pos * 3] >>> offset) & 0x01;
        int code = (int) (seq[pos * 3 + (offset >> 5) + 1] >>> shift) & 0x03;
        return code + ((4 - code) * nFlag);
    }

    @Override
    public void set(long index, long val) {
        // |N64 ... N0|B0 B1 ....  B31|B32 B33 ... B63|
        int pos = (int) (index >> 6);
        int offset = (int) (index & 0x3FL);
        int shift = (offset & 0x1F) << 1;

        seq[pos * 3] &= ~(0x01L << offset);
        seq[pos * 3] |= ((val & 0x04L) >>> 2) << offset;
        int bPos = pos * 3 + (offset >> 5) + 1;
        seq[bPos] &= ~(0xC000000000000000L >>> shift);
        seq[bPos] |= (val & 0x03) << (62 - shift);
    }

    @Override
    public long textSize() {
        return numBases;
    }

    /**
     * Create a reverse string of the this sequence. The sentinel is appended as
     * the last character of the resulting sequence.
     * 
     * @return Reverse sequence. The returned sequence is NOT a complement of
     *         the original sequence.
     */
    public ACGTSequence reverse() {
        ACGTSequence rev = new ACGTSequence(this.numBases);
        for (long i = 0; i < numBases; ++i) {
            rev.set(i, this.lookup(numBases - i - 1));
        }
        return rev;
    }

    /**
     * Create a complementary sequence (not reversed)
     * 
     * @return complementary sequence
     */
    public ACGTSequence complement() {
        ACGTSequence c = new ACGTSequence(this.numBases);

        int numBlocks = seq.length / 3;
        for (int i = 0; i < numBlocks; ++i) {
            c.seq[i * 3] = this.seq[i * 3];
            c.seq[i * 3 + 1] = ~(this.seq[i * 3 + 1]);
            c.seq[i * 3 + 2] = ~(this.seq[i * 3 + 2]);
        }
        return c;
    }

    public static ACGTSequence loadFrom(File f) throws IOException {
        DataInputStream d = new DataInputStream(new BufferedInputStream(new FileInputStream(f), 4 * 1024 * 1024));
        try {
            return ACGTSequence.loadFrom(d);
        }
        finally {
            d.close();
        }
    }

    public static ACGTSequence loadFrom(DataInputStream in) throws IOException {
        // The num bases must be always 2

        long numBases = in.readLong();
        int longArraySize = (int) ((numBases * 3 + LONG_BYTE_SIZE - 1) / LONG_BYTE_SIZE);
        long[] seq = new long[longArraySize];
        SnappyInputStream sin = new SnappyInputStream(in);
        int readBytes = sin.read(seq);
        return new ACGTSequence(seq, numBases);
    }

    public void saveTo(DataOutputStream out) throws IOException {
        out.writeLong(this.numBases);
        SnappyOutputStream sout = new SnappyOutputStream(out);
        sout.write(seq);
        sout.flush();
    }

    public void saveTo(File file) throws IOException {
        DataOutputStream d = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file)));
        try {
            saveTo(d);
        }
        finally {
            d.close();
        }
    }

    @Override
    public String toString() {
        StringBuilder b = new StringBuilder();
        for (int i = 0; i < numBases; ++i) {
            ACGT base = ACGT.decode((byte) lookup(i));
            b.append(base);
        }
        return b.toString();
    }

    /**
     * Count the number of the specified character in the range
     * 
     * @param code
     * @param start
     * @param end
     * @return
     */
    public long count(char code, long start, long end) {
        long count = 0;
        for (long i = start; i < end; ++i) {
            if ((char) lookup(i) == code)
                count++;
        }
        return count;
    }

    /**
     * Count the number of occurrence of the code within the specified range
     * 
     * @param code
     * @param start
     * @param end
     * @return
     */
    public long fastCount(IUPAC code, long start, long end) {
        long count = 0;
        long cursor = start;
        if (cursor < end && cursor % 2 != 0) {
            if (lookup(cursor) == code.bitFlag)
                count++;
            cursor++;
        }

        for (; cursor + 16 < end; cursor += 16) {
            int pos = (int) (cursor >>> 1);
            long v = 0;
            // Fill a long value from the byte array [pos, ... pos+8)
            for (int i = 0; i < 8; ++i) {
                v <<= 8;
                v |= seq[pos + i] & 0xFF;
            }

            long r = ~0L;
            for (int k = 0; k < 4; ++k) {
                r &= (((code.bitFlag & (0x08 >>> k)) == 0 ? ~v : v) << k) & 0x8888888888888888L;
            }
            count += countOneBit(r);
        }

        for (; cursor < end; cursor++) {
            if (lookup(cursor) == code.bitFlag)
                count++;
        }
        return count;
    }

    static int interleaveWith0(int v) {
        v = ((v & 0xFF00) << 8) | (v & 0x00FF);
        v = ((v << 4) | v) & 0x0F0F0F0F;
        v = ((v << 2) | v) & 0x33333333;
        v = ((v << 1) | v) & 0x55555555;
        return v;
    }

    static long interleave32With0(long v) {
        v = ((v & 0xFFFF0000L) << 16) | (v & 0x0000FFFFL);// 0000000000000000
        v = ((v << 8) | v) & 0x00FF00FF00FF00FFL; // 0000000011111111
        v = ((v << 4) | v) & 0x0F0F0F0F0F0F0F0FL; // 00001111
        v = ((v << 2) | v) & 0x3333333333333333L; // 0011
        v = ((v << 1) | v) & 0x5555555555555555L; // 0101
        return v;
    }

    /**
     * Count the number of 1s in the input. See also the Hacker's Delight:
     * http://hackers-delight.org.ua/038.htm
     * 
     * @param x
     * @return the number of 1-bit in the input x
     */
    public static long countOneBit(long x) {
        x = (x & 0x5555555555555555L) + ((x >>> 1) & 0x5555555555555555L);
        x = (x & 0x3333333333333333L) + ((x >>> 2) & 0x3333333333333333L);
        x = (x + (x >>> 4)) & 0x0F0F0F0F0F0F0F0FL;
        x = x + (x >>> 8);
        x = x + (x >>> 16);
        x = x + (x >>> 32);
        return x & 0x7FL;
    }

    @Override
    public long increment(long i, long val) {
        throw new UnsupportedOperationException("update");
    }

}
