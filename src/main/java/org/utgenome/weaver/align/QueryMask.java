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
// QueryMask.java
// Since: 2011/08/15
//
// $URL$ 
// $Author$
//--------------------------------------package org.utgenome.weaver.align.strategy;
package org.utgenome.weaver.align;

/**
 * A set of bit flags of ACGT characters in a query sequence
 * 
 * @author leo
 * 
 */
public class QueryMask
{
    private int         m;
    private BitVector[] patternMaskF;
    private BitVector[] patternMaskR; // reverse pattern

    public QueryMask(ACGTSequence query) {
        m = (int) query.textSize();
        patternMaskF = new BitVector[ACGT.exceptN.length];
        patternMaskR = new BitVector[ACGT.exceptN.length];
        for (int i = 0; i < patternMaskF.length; ++i) {
            patternMaskF[i] = new BitVector(m);
            patternMaskR[i] = new BitVector(m);
        }

        for (int i = 0; i < m; ++i) {
            ACGT ch = query.getACGT(i);
            if (ch == ACGT.N) {
                for (ACGT each : ACGT.exceptN) {
                    patternMaskF[each.code].set(i);
                    patternMaskR[each.code].set(m - i - 1);
                }
            }
            else {
                patternMaskF[ch.code].set(i);
                patternMaskR[ch.code].set(m - i - 1);
            }
        }
    }

    /**
     * Get bit flags of 64-bit range of the specified character in the query
     * sequence
     * 
     * @param ch
     * @param start
     * @return
     */
    public long getPatternMaskIn64bitForBidirectionalSearch(ACGT ch, int start) {
        if (start < 0) {
            long p = patternMaskF[ch.code].substring64(0, 64);
            return p << (-start);
        }
        long p = patternMaskF[ch.code].substring64(start, m);
        if (start + 64 >= m) {
            // combine forward and reverse pattern mask
            p |= (patternMaskR[ch.code].substring64(0, start + 64 - m)) << (m - start);
        }
        return p;
    }

    /**
     * Get bit flags of 64-bit range of the specified character in the query
     * sequence
     * 
     * @param ch
     * @param start
     * @return
     */
    public long getPatternMaskIn64bit(ACGT ch, int blockIndex, int w) {
        return patternMaskF[ch.code].substring64(blockIndex * w, (blockIndex + 1) * w);
    }

    @Override
    public String toString() {
        StringBuilder s = new StringBuilder();
        int i = 0;
        for (; i < ACGT.exceptN.length; ++i) {
            if (i > 0)
                s.append(", ");
            s.append(String.format("%s:%s", ACGT.exceptN[i], patternMaskF[i]));
        }
        return s.toString();
    }
}