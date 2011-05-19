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
// IUPACSequenceReader.java
// Since: 2011/02/10
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

import org.utgenome.gwt.utgb.client.bio.IUPAC;
import org.xerial.util.log.Logger;

/**
 * Reader of the IUPACSequence
 * 
 * @author leo
 * 
 */
public class IUPACSequence implements LSeq
{
    private static Logger _logger = Logger.getLogger(IUPACSequence.class);

    private byte[]        seq;
    private long          numBases;

    public IUPACSequence(String s) {
        this(s, false);
    }

    public IUPACSequence(String s, boolean appendSentinel) {
        this.numBases = s.length() + 1; // add 1 for sentinel
        long byteSize = (numBases / 2) + (numBases & 0x01);
        seq = new byte[(int) byteSize];
        for (int i = 0; i < seq.length; ++i)
            seq[i] = 0;

        int cursor = 0;
        for (; cursor < s.length(); ++cursor) {
            char c = Character.toUpperCase(s.charAt(cursor));
            setIUPAC(cursor, IUPAC.encode(c));
        }
        if (appendSentinel) {
            if (cursor % 2 == 0)
                setIUPAC(cursor++, IUPAC.N);
            // append a sentinal
            setIUPAC(cursor++, IUPAC.None);
        }

        numBases = cursor;
    }

    public IUPACSequence(long numBases) {
        this.numBases = numBases;
        long byteSize = (numBases / 2) + (numBases & 0x01);
        if (byteSize > Integer.MAX_VALUE)
            throw new IllegalArgumentException(String.format("iupac sequence cannot be larger than 4GB: %,d", numBases));

        this.seq = new byte[(int) byteSize];

        for (int i = 0; i < seq.length; ++i)
            seq[i] = 0;
    }

    private IUPACSequence(long numBases, byte[] seq) {
        this.numBases = numBases;
        this.seq = seq;
    }

    /**
     * Create a reverse string of the this sequence. The sentinel is appended as
     * the last character of the resulting sequence.
     * 
     * @return Reverse sequence. The returend sequence is not a complement of
     *         the original sequence.
     */
    public IUPACSequence reverse() {
        IUPACSequence rev = new IUPACSequence(numBases);
        long cursor = 0;
        long i = numBases - 1;
        boolean hasSentinel = false;
        if (getIUPAC(textSize() - 1) == IUPAC.None) {
            // Ignore the sentinel in the last character
            hasSentinel = true;
            i--;
        }

        // Reverse the sequence 
        for (; i >= 0; --i) {
            rev.setIUPAC(cursor++, getIUPAC(i));
        }
        if (hasSentinel) {
            // Append a sentinel
            rev.setIUPAC(cursor, IUPAC.None);
        }
        return rev;
    }

    /**
     * Create a complementary sequence (not reversed)
     * 
     * @return complementary sequence
     */
    public IUPACSequence complement() {
        try {
            ByteArrayOutputStream out = new ByteArrayOutputStream(this.seq.length);
            IUPACSequenceWriter buf = new IUPACSequenceWriter(out);
            for (long i = 0; i < textSize(); ++i) {
                buf.append(getIUPAC(i).complement());
            }
            buf.close();
            return new IUPACSequence(textSize(), out.toByteArray());
        }
        catch (IOException e) {
            throw new IllegalStateException("Unexpected IOException: " + e.getMessage());
        }
    }

    public static IUPACSequence loadFrom(File f) throws IOException {
        long numBases = f.length() * 2;
        DataInputStream d = new DataInputStream(new BufferedInputStream(new FileInputStream(f), 4 * 1024 * 1024));
        try {
            return IUPACSequence.loadFrom(d, numBases);
        }
        finally {
            d.close();
        }
    }

    public static IUPACSequence loadFrom(DataInputStream in, long numBases) throws IOException {
        // The num bases must be always 2
        int byteSize = (int) (numBases / 2);
        byte[] seq = new byte[byteSize];
        int readBytes = in.read(seq, 0, byteSize);

        return new IUPACSequence(numBases, seq);
    }

    public void saveTo(DataOutputStream out) throws IOException {
        out.write(seq, 0, seq.length);
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

    public IUPAC getIUPAC(long index) {
        byte code = (byte) ((seq[(int) (index >>> 1)] >>> (1 - (index & 1)) * 4) & 0x0F);
        return IUPAC.decode(code);
    }

    public void setIUPAC(long index, IUPAC val) {
        int pos = (int) (index / 2);
        int offset = (int) (index % 2);
        byte code = (byte) val.bitFlag;

        byte mask = (byte) ~(0xF0 >>> (offset * 4));
        seq[pos] &= mask;
        seq[pos] |= (byte) (code << (1 - offset) * 4);
    }

    @Override
    public String toString() {
        StringBuilder b = new StringBuilder();
        for (int i = 0; i < numBases; ++i) {
            IUPAC base = getIUPAC(i);
            b.append(base == IUPAC.None ? "$" : base.name());
        }
        return b.toString();
    }

    public String toACGTString() {
        StringBuilder b = new StringBuilder();
        for (int i = 0; i < numBases; ++i) {
            IUPAC base = getIUPAC(i);
            switch (base) {
            case A:
            case C:
            case G:
            case T:
                break;
            default:
                base = IUPAC.N;
            }
            b.append(base == IUPAC.None ? "$" : base.name());
        }
        return b.toString();
    }

    public long count(IUPAC code, long start, long end) {
        long count = 0;
        for (long i = start; i < end; ++i) {
            if (lookup(i) == code.bitFlag)
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
    public long lookup(long index) {
        int pos = (int) (index >>> 1);
        int shift = 4 * (1 - (int) (index & 1));
        return ((seq[pos] >>> shift) & 0x0F);
    }

    @Override
    public long textSize() {
        return numBases;
    }

    @Override
    public void set(long i, long val) {
        throw new UnsupportedOperationException("set");
    }

    @Override
    public long update(long i, long val) {
        throw new UnsupportedOperationException("update");
    }

}
