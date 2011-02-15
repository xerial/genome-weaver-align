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
// Since: 2011/02/15
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

/**
 * Bit-vector supporting rank/select operations
 * 
 * @author leo
 * 
 */
public class BitVector
{
    private static final int NBITS_IN_LONG            = 64;
    private static final int NUM_BLOCKS_IN_RANK_TABLE = 4;
    private LLongArray       block;
    private LLongArray       rankTable;
    private long             size;

    public BitVector(long size) {
        this.size = size;
        long numBlocks = (size + NBITS_IN_LONG - 1) / NBITS_IN_LONG;
        block = new LLongArray(numBlocks);
        clear();
    }

    public long size() {
        return size;
    }

    public void clear() {
        block.fill(0L);
        rankTable = null;
    }

    public void prepareRankTable() {
        long tableSize = (block.size() + NUM_BLOCKS_IN_RANK_TABLE - 1) / NUM_BLOCKS_IN_RANK_TABLE + 1;
        rankTable = new LLongArray(tableSize);
        rankTable.fill(0L);
        long count = 0;
        for (long i = 0; i < block.size(); ++i) {
            if ((i % NUM_BLOCKS_IN_RANK_TABLE) == 0) {
                rankTable.set(i / NUM_BLOCKS_IN_RANK_TABLE, count);
            }
            count += popCount(block.get(i));
        }
        rankTable.set(rankTable.size() - 1, count);
    }

    public boolean get(long index) {
        long blockPos = index / NBITS_IN_LONG;
        if (blockPos > block.size())
            return false;

        long offset = index % NBITS_IN_LONG;
        long mask = 1L << offset;
        return (block.get(blockPos) & mask) != 0;
    }

    public void set(long index) {
        long blockPos = index / NBITS_IN_LONG;
        long offset = index % NBITS_IN_LONG;
        block.setOR(blockPos, 1L << offset);
    }

    public void reset(long index) {
        long blockPos = index / NBITS_IN_LONG;
        long offset = index % NBITS_IN_LONG;
        long mask = 1L << offset;
        block.setAND(blockPos, ~mask);
    }

    public long rank(boolean c, long index) {
        if (index > size)
            return 0;
        if (c)
            return rankOne(index);
        else
            return index - rankOne(index);
    }

    private long rankOne(long index) {
        if (rankTable == null) {
            prepareRankTable();
        }

        long blockPos = index / NBITS_IN_LONG;
        long tablePos = blockPos / NUM_BLOCKS_IN_RANK_TABLE;

        long rank = rankTable.get(tablePos);
        for (long i = tablePos * NUM_BLOCKS_IN_RANK_TABLE; i < blockPos; ++i) {
            rank += popCount(block.get(i));
        }
        long offset = index % NBITS_IN_LONG;
        if (offset != 0) {
            long val = block.get(blockPos);
            long mask = (1L << offset) - 1;
            rank += popCount(val & mask);
        }
        return rank;
    }

    /**
     * Count the number of 1s in the input
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
        return x & 0x7F;
    }

    @Override
    public String toString() {
        StringBuilder buf = new StringBuilder();
        for (int i = 0; i < size; i++)
            buf.append(get(i) ? "1" : "0");
        return buf.toString();
    }

}
