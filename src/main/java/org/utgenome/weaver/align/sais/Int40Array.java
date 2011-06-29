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
// UInt40Array.java
// Since: 2011/05/19
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.sais;

import java.util.Iterator;

import org.utgenome.weaver.align.LSeq;

/**
 * ArrayList for int40 values ([-(2^39-1), 2^39-1])
 * 
 * TODO implementation
 * 
 * @author leo
 * 
 */
public class Int40Array implements LSeq, Iterable<Long>
{
    private static final int  B               = 30;                // bit length 
    private static final int  CONTAINER_SIZE  = 1 << B;
    private static final long OFFSET_MASK     = CONTAINER_SIZE - 1;
    private static final int  LONG_BYTE_SIZE  = 8;
    private static final int  INT40_BYTE_SIZE = 5;
    private final long        size;

    private long[]            rawArray;

    public Int40Array(long size) {
        this.size = size;
        long byteSize = size * INT40_BYTE_SIZE;

        // 8 * 2GB = 16GB 
        // int40   16GB / 5 = 3.2G entries
        // uint32  16GB / 4 = 4G entries
        long rawArraySize = (byteSize + LONG_BYTE_SIZE - 1) / LONG_BYTE_SIZE;
        if (rawArraySize > Integer.MAX_VALUE)
            throw new IllegalArgumentException(String.format(
                    "Cannot create int40 array more than %,d entries: size=%,d", Integer.MAX_VALUE, size));
        rawArray = new long[(int) rawArraySize];
    }

    public long lookup(long index) {

        final long bytePos = index * INT40_BYTE_SIZE;
        final int pos = (int) (bytePos / LONG_BYTE_SIZE);
        final int offset = (int) (bytePos % LONG_BYTE_SIZE);

        final int maskEnd = offset + INT40_BYTE_SIZE;
        // byte index: 01234567
        // byte mask:  0XXXXX00
        //   | <- offset = 1, shiftLen = 8 - offset - 5= 2 
        //
        // byte index: 01234567
        // byte mask:  00000XXXXX00
        //   | <- offset = 5, shiftLen = 8 - offset - 5 = -2

        long v = rawArray[pos] & (0xFFFFFFFFFF000000L >>> (offset * 8));
        int shiftLen = LONG_BYTE_SIZE - maskEnd;
        if (shiftLen >= 0) {
            v >>>= shiftLen * 8;
        }
        else {
            v <<= -shiftLen * 8;
        }

        // byte index: 0123456701234567
        // byte mask:  000000XXXXX00000
        //                   | <- offset = 6, maskEnd = 6 + 5 = 11, nextMaskEnd = 11 % 8 = 3;
        //                     000XXXXX   shift len = 5 
        if (offset >= 4) {
            int nextMaskEnd = maskEnd % LONG_BYTE_SIZE;
            long v2 = rawArray[pos + 1] & (0xFFFFFFFFFF000000L << ((INT40_BYTE_SIZE - nextMaskEnd) * 8));
            v2 >>>= (LONG_BYTE_SIZE - nextMaskEnd) * 8;
            v |= v2;
        }

        if ((v & 0x8000000000L) == 0)
            return v;
        else
            return v - 0xFFFFFFFFFFL - 1;
    }

    public void set(long index, long value) {
        final long bytePos = index * INT40_BYTE_SIZE;
        final int pos = (int) (bytePos / LONG_BYTE_SIZE);
        final int offset = (int) (bytePos % LONG_BYTE_SIZE);

        // byte index: 01234567
        // value:      000XXXXX
        //
        // offset:          | <- offset 5, shiftLen = 2
        // shift       00000XXXXX
        // 
        // offset:      | <- offset 1, shiftLen = -2
        // shift       0XXXXX00
        long v = value & 0xFFFFFFFFFFL;

        long mask = 0xFFFFFFFFFF000000L >>> (offset * 8);
        rawArray[pos] &= ~mask;
        int shiftLen = offset - 3;
        rawArray[pos] |= shiftLen >= 0 ? v >>> (shiftLen * 8) : v << -shiftLen * 8;

        // byte index: 0123456701234567
        // value:      000XXXXX
        //
        // offset:          | <- offset 5 
        // shift       00000XXX
        // shift       XX000000
        if (offset >= 4) {
            int shiftLenOfSecondHalf = LONG_BYTE_SIZE - shiftLen;
            rawArray[pos + 1] &= ~(0xFFFFFFFFFFL << (shiftLenOfSecondHalf * 8));
            rawArray[pos + 1] |= v << (shiftLenOfSecondHalf * 8);
        }
    }

    public long textSize() {
        return size;
    }

    @Override
    public long increment(long index, long val) {
        long next = lookup(index) + val;
        set(index, next);
        return next;
    }

    @Override
    public String toString() {
        StringBuilder b = new StringBuilder();
        int i = 0;
        b.append("[");
        for (long each : this) {
            if (i++ > 0)
                b.append(", ");
            b.append(each);
        }
        b.append("]");
        return b.toString();
    }

    public long[] toArray() {
        long[] r = new long[(int) textSize()];
        for (long i = 0; i < textSize(); ++i) {
            r[(int) i] = lookup(i);
        }
        return r;
    }

    @Override
    public Iterator<Long> iterator() {
        return new Iterator<Long>() {

            private long cursor = 0;

            @Override
            public boolean hasNext() {
                return cursor < size;
            }

            @Override
            public Long next() {
                return lookup(cursor++);
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException("remove");
            }
        };
    }

}
