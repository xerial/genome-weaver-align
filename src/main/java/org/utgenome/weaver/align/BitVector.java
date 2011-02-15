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
    private static final int B          = 64; // the number of bits in a block
    private static final int M          = 4; // block width to precompute the ranks  
    private LLongArray       block;
    private LLongArray       rankTable;
    private long             size;
    private long             numberOf1s = 0;

    public BitVector(long size) {
        this.size = size;
        long numBlocks = (size + B - 1) / B;
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
        long tableSize = (block.size() + M - 1) / M + 1;
        rankTable = new LLongArray(tableSize);
        rankTable.fill(0L);
        long count = 0;
        for (long i = 0; i < block.size(); ++i) {
            if ((i % M) == 0) {
                rankTable.set(i / M, count);
            }
            count += popCount(block.get(i));
        }
        rankTable.set(rankTable.size() - 1, count);
        numberOf1s = count;
    }

    public boolean get(long index) {
        long blockPos = index / B;
        if (blockPos > block.size())
            return false;

        long offset = index % B;
        long mask = 1L << offset;
        return (block.get(blockPos) & mask) != 0;
    }

    public void setBit(boolean c, long index) {
        if (c)
            set(index);
        else
            reset(index);
    }

    public void set(long index) {
        long blockPos = index / B;
        long offset = index % B;
        block.setOR(blockPos, 1L << offset);
    }

    public void reset(long index) {
        long blockPos = index / B;
        long offset = index % B;
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
        if (rankTable == null)
            prepareRankTable();

        long blockPos = index / B;
        long tablePos = blockPos / M;

        long rank = rankTable.get(tablePos);
        for (long i = tablePos * M; i < blockPos; ++i) {
            rank += popCount(block.get(i));
        }
        long offset = index % B;
        if (offset != 0) {
            long val = block.get(blockPos);
            long mask = (1L << offset) - 1;
            rank += popCount(val & mask);
        }
        return rank;
    }

    public long select(boolean c, long rank) {
        if (rankTable == null)
            prepareRankTable();

        if (c) {
            if (rank > numberOf1s)
                return -1;
        }
        else {
            if (rank > size - numberOf1s)
                return -1;
        }

        RankAndPos rp = selectOutBlock(c, rank);
        long val = c ? block.get(rp.pos) : ~block.get(rp.pos);
        return rp.pos * B + selectInBlock(val, rp.rank);
    }

    private long getBitNum(long oneNum, long num, boolean bit) {
        if (bit)
            return oneNum;
        else
            return num - oneNum;
    }

    private static class RankAndPos
    {
        public final long rank;
        public final long pos;

        public RankAndPos(long rank, long pos) {
            this.rank = rank;
            this.pos = pos;
        }
    }

    private RankAndPos selectOutBlock(boolean c, long rank) {
        // binary search over tables
        long left = 0;
        long right = rankTable.size();
        while (left < right) {
            long mid = (left + right) / 2;
            long length = B * M * mid;
            if (getBitNum(rankTable.get(mid), length, c) < rank) {
                left = mid + 1;
            }
            else {
                right = mid;
            }
        }

        long table_ind = (left != 0) ? left - 1 : 0;
        long block_pos = table_ind * M;
        rank -= getBitNum(rankTable.get(table_ind), block_pos * B, c);

        // sequential search over blocks
        for (; block_pos < block.size(); ++block_pos) {
            long rank_next = getBitNum(popCount(block.get(block_pos)), B, c);
            if (rank <= rank_next) {
                break;
            }
            rank -= rank_next;
        }

        return new RankAndPos(rank, block_pos);
    }

    private long selectInBlock(long x, long rank) {
        long x1 = x - ((x & 0xAAAAAAAAAAAAAAAAL) >>> 1);
        long x2 = (x1 & 0x3333333333333333L) + ((x1 >>> 2) & 0x3333333333333333L);
        long x3 = (x2 + (x2 >>> 4)) & 0x0F0F0F0F0F0F0F0FL;

        long pos = 0;
        for (;; pos += 8) {
            long rank_next = (x3 >>> pos) & 0xFFL;
            if (rank <= rank_next)
                break;
            rank -= rank_next;
        }

        long v2 = (x2 >>> pos) & 0xFL;
        if (rank > v2) {
            rank -= v2;
            pos += 4;
        }

        long v1 = (x1 >> pos) & 0x3L;
        if (rank > v1) {
            rank -= v1;
            pos += 2;
        }

        long v0 = (x >> pos) & 0x1L;
        if (v0 < rank) {
            rank -= v0;
            pos += 1;
        }

        return pos;

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

}
