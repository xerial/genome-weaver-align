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
    private ACGT              leadingBase    = ACGT.T;
    private long[]            seq;
    private long              numBases;
    private long              capacity;

    public SOLiDColorSequence() {
        ensureArrayCapacity(10);
        numBases = 0;
    }

    public SOLiDColorSequence(String s) {
        this(s.length() - 1);

        if (s.length() <= 0) {
            this.leadingBase = ACGT.N;
        }
        else {
            this.leadingBase = ACGT.encode(s.charAt(0));
            for (int i = 1; i < s.length(); ++i) {
                SOLiDColor code = SOLiDColor.encode(s.charAt(i));
                set(i - 1, code);
            }
        }
    }

    public SOLiDColorSequence(ACGTSequence acgt) {
        this(acgt.textSize());

        ACGT prev = leadingBase;
        for (long i = 0; i < acgt.textSize(); ++i) {
            ACGT next = acgt.getACGT(i);
            SOLiDColor code = SOLiDColor.encode(prev, next);
            set(i, code);
            prev = next;
        }
    }

    public SOLiDColorSequence(long numBases) {
        this.numBases = numBases;

        ensureArrayCapacity(numBases);
    }

    public ACGTSequence toACGTSequence() {
        ACGTSequence acgt = new ACGTSequence(numBases);
        ACGT current = leadingBase;
        for (long i = 0; i < numBases; ++i) {
            current = SOLiDColor.decode(current, getColor(i));
            acgt.set(i, current);
        }
        return acgt;
    }

    public String toColorString() {
        StringBuilder s = new StringBuilder((int) numBases);
        for (int i = 0; i < numBases; ++i) {
            s.append(getColor(i));
        }
        return s.toString();
    }

    public SOLiDColorSequence reverseComplement() {
        ACGTSequence rc = this.toACGTSequence().reverseComplement();
        SOLiDColorSequence rev = new SOLiDColorSequence(this.numBases);
        if (rc.textSize() > 0 && rc.getACGT(0) != ACGT.N) {
            // Use the default leading base
            rev.set(0, SOLiDColor.encode(rev.leadingBase, rc.getACGT(0)));
        }
        else {
            // Use N when leading base is unknown
            rev.leadingBase = ACGT.N;
        }
        for (int i = 0; i < this.numBases; ++i) {
            rev.set(i + 1, getColor(numBases - i - 1));
        }
        return rev;
    }

    @Override
    public String toString() {
        return String.format("%s%s", leadingBase, toColorString());
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

    public void set(long index, SOLiDColor color) {
        set(index, color.code);
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
