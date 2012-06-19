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
 * ArrayList for int40 values. This array can hold 3.2G entries.
 * 
 * 
 * @author leo
 * 
 */
public class Int40Array implements LSeq, Iterable<Long>
{
    public static final long MAX_VALUE       = 0x7FFFFFFFFFL;
    public static final long MIN_VALUE       = -0x7FFFFFFFFFL;

    public static final long MAX_INDEX       = 16L * 1024L * 1024L * 1024L / 5L;

    private static final int LONG_BYTE_SIZE  = 8;
    private static final int INT40_BYTE_SIZE = 5;
    private final long       size;

    private long[]           rawArray;

    public Int40Array(long size) {
        this.size = size;

        long numBlocks = (size * 5 + 5 - 1) / 5;
        // 8 * 2GB = 16GB 
        // int40   16GB / 5 = 3.2G entries
        // uint32  16GB / 4 = 4G entries
        long rawArraySize = numBlocks * 5;
        if (rawArraySize > MAX_INDEX)
            throw new IllegalArgumentException(String.format(
                    "Cannot create int40 array more than %,d entries: size=%,d", MAX_INDEX, size));
        rawArray = new long[(int) rawArraySize];
    }

    public long lookup(long index) {
        int bPos = (int) index >> 2;
        int bOffset = (int) index & 0x03;
        long L = rawArray[bPos * 5 + bOffset + 1] >>> ((1 - (bOffset & 0x01)) * 32);
        long H = ((rawArray[bPos * 5] >>> (bOffset * 8)) & 0xFF) << 32;
        long v = H | L;
        return v - ((v & 0x8000000000L) << 1);
    }

    public void set(long index, long value) {
        int bPos = (int) index >> 2;
        int bOffset = (int) index & 0x03;
        long L = value & 0xFFFFFFFFL;
        long H = value >>> 32;
        rawArray[bPos * 5 + bOffset + 1] &= 0xFFFFFFFFL >>> ((index & 0x01) * 32);
        rawArray[bPos * 5 + bOffset + 1] |= L << ((1 - index & 0x01) * 32);
        rawArray[bPos * 5] &= ~(0xff << (bOffset * 8));
        rawArray[bPos * 5] |= H << (bOffset * 8);
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
