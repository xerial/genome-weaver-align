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
// BitVector.java
// Since: 2011/07/20
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;

/**
 * Bit vector of 0/1 binary string
 * 
 * @author leo
 * 
 */
public class BitVector
{
    private static final int B = 64; // the number of bits in a block
    private final long[]     block;
    private final long       size;

    public BitVector(long size) {
        this.size = size;
        long blockSize = (size + B - 1) / B;
        if (blockSize > Integer.MAX_VALUE)
            throw new IllegalArgumentException("cannot create bit vector with size > " + Integer.MAX_VALUE
                    + ". size = " + size);

        this.block = new long[(int) blockSize];
        clear();
    }

    private BitVector(long size, long[] block) {
        this.size = size;
        this.block = block;
    }

    public long size() {
        return size;
    }

    public void clear() {
        for (int i = 0; i < block.length; ++i)
            block[i] = 0L;
    }

    public boolean get(long index) {
        long blockPos = index / B;
        if (blockPos > block.length)
            return false;

        long offset = index % B;
        long mask = 1L << offset;
        return (block[(int) blockPos] & mask) != 0;
    }

    public int getInt(long index) {
        long blockPos = index / B;
        if (blockPos > block.length)
            return 0;

        long offset = index % B;
        long mask = 1L << offset;
        return (int) ((block[(int) blockPos] & mask) & 0x1L);
    }

    public void setBit(boolean c, long index) {
        if (c)
            set(index);
        else
            reset(index);
    }

    /**
     * Set the bit at the index. Note that this operation invalidates the
     * current rank table, so you must refresh the rank table by calling
     * {@link #refreshRankTable()}
     * 
     * @param index
     */
    public void set(long index) {
        long blockPos = index / B;
        long offset = index % B;
        block[(int) blockPos] |= 1L << offset;
    }

    public void reset(long index) {
        long blockPos = index / B;
        long offset = index % B;
        long mask = 1L << offset;
        block[(int) blockPos] &= ~mask;
    }

    /**
     * Count the number of 1s in the input. See also the Hacker's Delight:
     * http://hackers-delight.org.ua/038.htm
     * 
     * @param x
     * @return the number of 1-bit in the input x
     */
    public static long popCount(long x) {
        x = (x & 0x5555555555555555L) + ((x >>> 1) & 0x5555555555555555L);
        x = (x & 0x3333333333333333L) + ((x >>> 2) & 0x3333333333333333L);
        x = (x + (x >>> 4)) & 0x0F0F0F0F0F0F0F0FL;
        x = x + (x >>> 8);
        x = x + (x >>> 16);
        x = x + (x >>> 32);
        return x & 0x7FL;
    }

    @Override
    public String toString() {
        StringBuilder buf = new StringBuilder();
        for (int i = 0; i < size; i++)
            buf.append(get(i) ? "1" : "0");
        return buf.toString();
    }

    public DataOutputStream saveTo(DataOutputStream out) throws IOException {
        out.writeLong(size);
        out.writeInt(block.length);
        for (int i = 0; i < block.length; ++i)
            out.writeLong(block[i]);
        return out;
    }

    public static BitVector loadFrom(DataInputStream in) throws IOException {
        long size = in.readLong();
        int blockSize = in.readInt();
        long[] block = new long[blockSize];
        for (int i = 0; i < blockSize; ++i)
            block[i] = in.readLong();
        BitVector v = new BitVector(size, block);
        return v;
    }

}
