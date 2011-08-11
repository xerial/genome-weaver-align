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

public class StaircaseFilter
{
    private final int   m;
    private final int   k;

    private byte[]      chunkStart;
    private byte[]      chunkLen;
    private BitVector[] stairMask;

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

        stairMask = new BitVector[k + 1];
        for (int i = 0; i < stairMask.length; ++i) {
            stairMask[i] = new BitVector(m)._not()._lshift(chunkStart[i]);
        }
    }

    public BitVector getStaircaseMask(int k) {
        return stairMask[k];
    }

    public long getStairCaseMask64bit(int k, int offset) {
        if (offset >= 0)
            return stairMask[k].substring64(offset, offset + 64);
        else {
            return stairMask[k].substring64(0, 64) << (-offset);
        }
    }

    @Override
    public String toString() {
        return StringUtil.join(stairMask, ", ");
    }
}
