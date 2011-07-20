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
// SOLiDColorSequence.java
// Since: 2011/07/20
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.util.Arrays;

/**
 * SOLiD color sequence with N, using 3-bit encoding
 * 
 * @author leo
 * 
 */
public class SOLiDColorSequence implements LSeq
{
    private static final int  LONG_BYTE_SIZE = 8;
    private static final long MAX_SIZE       = (2 * 1024 * 1024 * 1024 * 8 * 8) / 3;

    // 2G (max of Java array size) * 8 (long byte size) * 8 / 3 (bit) =  42G  (42G characters)
    private ACGT              leadingBase;
    private long[]            seq;
    private long              numBases;
    private long              capacity;

    public SOLiDColorSequence() {
        ensureArrayCapacity(10);
        numBases = 0;
    }

    public SOLiDColorSequence(String s) {
        this(s.length());

        for (int i = 0; i < s.length(); ++i) {
            byte code = ACGT.to3bitCode(s.charAt(i));
            set(i, code);
        }
    }

    public SOLiDColorSequence(long numBases) {
        this.numBases = numBases;

        ensureArrayCapacity(numBases);
    }

    private static int minArraySize(long numBases) {
        long bitSize = numBases * 3L;
        long blockBitSize = LONG_BYTE_SIZE * 3L * 8L;
        long arraySize = ((bitSize + blockBitSize - 1L) / blockBitSize) * 3L;
        if (arraySize > Integer.MAX_VALUE) {
            throw new IllegalArgumentException(String.format("Cannot create ACGTSequece more than %,d size: %,d",
                    MAX_SIZE, numBases));
        }
        return (int) arraySize;
    }

    private void ensureArrayCapacity(long newCapacity) {
        if (seq != null && newCapacity < (seq.length * 64L / 3L)) {
            return;
        }
        long arraySize = minArraySize(newCapacity);

        if (seq == null) {
            seq = new long[(int) arraySize];
        }
        else {
            seq = Arrays.copyOf(seq, (int) arraySize);
        }
        this.capacity = (seq.length / 3L) * 64L;
    }

    private SOLiDColorSequence(long[] rawSeq, long numBases) {
        this.seq = rawSeq;
        this.numBases = numBases;
    }

    public SOLiDColor getColor(long index) {
        return SOLiDColor.decode((int) lookup(index));
    }

    @Override
    public long lookup(long index) {
        int pos = (int) (index >> 6);
        int offset = (int) (index & 0x03FL);
        int shift = 62 - ((int) (index & 0x1FL) << 1);

        long nFlag = seq[pos * 3] & (1L << (63 - offset));
        int code = (int) (seq[pos * 3 + (offset >> 5) + 1] >>> shift) & 0x03;
        return nFlag == 0 ? code : 4;
    }

    @Override
    public void set(long index, long val) {
        // |N0 ... N63|B0 B1 ....  B31|B32 B33 ... B63|
        int pos = (int) (index >> 6);
        int offset = (int) (index & 0x3FL);
        int shift = (offset & 0x1F) << 1;

        seq[pos * 3] &= ~(1L << (63 - offset));
        seq[pos * 3] |= ((val >>> 2) & 0x01L) << (63 - offset);
        int bPos = pos * 3 + (offset >> 5) + 1;
        seq[bPos] &= ~(0xC000000000000000L >>> shift);
        seq[bPos] |= (val & 0x03) << (62 - shift);
    }

    public void append(SOLiDColor base) {
        long index = this.numBases++;
        if (index >= this.capacity) {
            long newCapacity = (index * 3L / 2L) + 64L;
            ensureArrayCapacity(newCapacity);
        }
        set(index, base.code);
    }

    @Override
    public long textSize() {
        return numBases;
    }

    @Override
    public long increment(long index, long value) {
        throw new UnsupportedOperationException("update");
    }

}
