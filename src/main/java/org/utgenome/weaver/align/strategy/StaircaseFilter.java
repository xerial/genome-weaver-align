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
// StaircaseFilter.java
// Since: 2011/08/09
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.strategy;

import org.utgenome.weaver.align.BitVector;
import org.xerial.util.StringUtil;

/**
 * Staircase filter is a mask over NFA for read alignment
 * 
 * @author leo
 * 
 */
public class StaircaseFilter
{
    private final int   m;
    private final int   k;

    private byte[]      chunkStart;
    private byte[]      chunkLen;
    private BitVector[] staircaseMask;

    public StaircaseFilter(int m, int k) {
        this.m = m;
        this.k = k;
        int lastChunkSize = (m - k >= 6) ? m * 2 / (k + 2) : m - k;

        // prefix length of each chunk of the filter
        chunkStart = new byte[k + 2];
        byte rest = (byte) (m - lastChunkSize);
        if (k == 0)
            chunkStart[0] = rest;
        else {
            for (int i = 0; i <= k; ++i) {
                chunkStart[i] = (byte) (rest * i / k);
            }
        }
        chunkStart[k + 1] = (byte) m;

        chunkLen = new byte[k + 1];
        for (int i = 0; i <= k; ++i) {
            chunkLen[i] = (byte) (chunkStart[i + 1] - chunkStart[i]);
        }

        staircaseMask = new BitVector[k + 1];
        for (int i = 0; i < staircaseMask.length; ++i) {
            staircaseMask[i] = new BitVector(m)._not()._lshift(chunkStart[i]);
        }
    }

    public int getNumChunks() {
        return chunkLen.length;
    }

    public int getChunkStart(int chunkIndex) {
        return chunkStart[chunkIndex];
    }

    public int getChunkSize(int chunkIndex) {
        return chunkLen[chunkIndex];
    }

    public BitVector getStaircaseMask(int k) {
        return staircaseMask[k];
    }

    public long getStairCaseMask64bit(int k, int offset) {
        if (k >= staircaseMask.length)
            return 0L;
        int w = m - offset;
        if (offset >= 0) {
            long base = ~0L << (m - offset);
            return base | staircaseMask[k].substring64(offset, offset + 64);
        }
        else {
            return staircaseMask[k].substring64(0, 64) << (-offset);
        }
    }

    @Override
    public String toString() {
        return StringUtil.join(staircaseMask, ", ");
    }
}
