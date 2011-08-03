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
// BitParallelSmithWaterman.java
// Since: 2011/07/28
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.util.Arrays;

import org.xerial.util.log.Logger;

/**
 * Bit-parallel algorithm for Smith-Waterman alignment
 * 
 * <h3>References</h3>
 * 
 * <ul>
 * <li>Myers, G. (1999). A fast bit-vector algorithm for approximate string
 * matching based on dynamic programming. Journal of the ACM (JACM).</li>
 * <li>H. Hyyr&ouml; and G. Navarro. Faster Bit-parallel Approximate String
 * Matching In Proc. 13th Combinatorial Pattern Matching (CPM'2002), LNCS 2373
 * http://sunsite.dcc.uchile.cl/ftp/users/gnavarro/hn02.ps.gz</li>
 * </ul>
 * 
 * @author leo
 * 
 */
public class BitParallelSmithWaterman
{
    private static Logger _logger = Logger.getLogger(BitParallelSmithWaterman.class);

    public static void align64(ACGTSequence ref, ACGTSequence query, int numAllowedDiff) {
        // Preprocessing
        long[] pm = new long[ACGT.values().length];
        int m = Math.min(64, (int) query.textSize());
        for (int i = 0; i < m; ++i) {
            pm[query.getACGT(i).code] |= 1L << i;

        }
        long vp = ~0L;
        long vn = 0L;
        int diff = m;

        if (_logger.isDebugEnabled()) {
            for (ACGT ch : ACGT.exceptN) {
                _logger.debug("peq[%s]:%s", ch, toBinary(pm[ch.code], m));
            }
        }
        // Searching
        for (int j = 0; j < ref.textSize(); ++j) {
            long x = pm[ref.getACGT(j).code];
            long d0 = ((vp + (x & vp)) ^ vp) | x | vn;
            long hp = (vn | ~(vp | d0));
            long hn = vp & d0;
            x = (hp << 1);
            vn = x & d0;
            vp = (hn << 1) | ~(x | d0);
            // diff represents the last row (C[m, j]) of the DP matrix
            diff += (int) ((hp >>> m) & 1L);
            diff -= (int) ((hn >>> m) & 1L);
            if (_logger.isDebugEnabled()) {
                _logger.debug("[%s] j:%2d, diff:%2d %1s hp:%s, hn:%s, vp:%s, vn:%s, d0:%s", ref.getACGT(j), j, diff,
                        diff <= numAllowedDiff ? "*" : "", toBinary(hp, m), toBinary(hn, m), toBinary(vp, m),
                        toBinary(vn, m), toBinary(d0, m));
            }

        }
    }

    public static String toBinary(long v, int m) {
        StringBuilder s = new StringBuilder();
        for (int i = 0; i <= m; ++i) {
            s.append((v & (1L << i)) > 0 ? "1" : "0");
        }
        return s.toString();
    }

    public static class AlignBlocks
    {

        private static final int w = 64;                  // word size
        private static final int S = ACGT.values().length; // alphabet size

        private final int        k;
        private final int        m;
        private final int        bMax;
        private long[]           vp_in;
        private long[]           vp;
        private long[]           vn;
        private long[][]         peq;                     // [# of block][A, C, G, T]
        private int[]            D;                       // D[block]

        public AlignBlocks(ACGTSequence query, int k) {
            this.m = (int) query.textSize();
            this.k = k;
            bMax = (int) Math.ceil((double) m / w);
            vp = new long[bMax];
            vn = new long[bMax];
            vp_in = new long[bMax];
            peq = new long[bMax][S];

            // Peq bit-vector holds flags of the character occurrence positions in the query
            for (int i = 0; i < m; ++i) {
                peq[i / w][query.getACGT(i).code] |= 1L << (i % w);
            }
            // Fill the flanking region with 1s
            {
                int f = m % w;
                long mask = ~0L >>> f << f;
                for (int i = 0; i < S; ++i)
                    peq[bMax - 1][i] |= mask;
            }
        }

        public void clear() {
            Arrays.fill(vp, 0L);
            Arrays.fill(vn, 0L);
            Arrays.fill(vp_in, 0L);
            for (int i = 0; i < bMax; ++i)
                for (int j = 0; j < S; ++j)
                    peq[i][j] = 0L;
        }

        public void align(ACGTSequence ref, ACGTSequence query) {
            int b = (int) Math.ceil((double) k / w);
            for (int r = 0; r < b; ++r) {
                vp[r] = ~0L; // all 1s
                vn[r] = 0L; // all 0s
                D[r] = b * w;
            }

            final int N = (int) ref.textSize();
            for (int j = 0; j < N; ++j) {
                ACGT ch = ref.getACGT(j);
                int carry = 0;
                for (int r = 0; r < b; ++r) {
                    D[b] += (carry = alignBlock(ch, r, carry));
                }

                if (D[b] - carry <= k && b < bMax && (((peq[b + 1][ch.code] & 1L) != 0) | carry < 0)) {
                    b++;
                    vp[b] = vp_in[b];
                    vn[b] = 0L;

                    D[b] = D[b - 1] + w - carry + alignBlock(ch, b, carry);
                }
                else {
                    while (D[b] >= k + w) {
                        --b;
                    }
                }

            }

        }

        private int alignBlock(ACGT ch, int r, int hin) {
            long vp = this.vp[r];
            long vn = this.vn[r];
            long x = this.peq[r][ch.code];
            if (hin < 0)
                x |= 1L;
            long d0 = ((vp + (x & vp)) ^ vp) | x | vn;
            long hn = vp & d0;
            long hp = (vn | ~(vp | d0));

            int hout = 0;
            hout += (int) ((hp >>> w) & 1L);
            hout -= (int) ((hn >>> w) & 1L);

            long ph = (hp << 1);
            long mh = (hn << 1);
            if (hin < 0)
                mh |= 1;
            else if (hin > 0)
                ph |= 1;

            this.vp[r] = (mh << 1) | ~(ph | d0);
            this.vn[r] = ph & d0;

            if (_logger.isDebugEnabled()) {
                _logger.debug("block:%d, hin:%2d, hout:%2d %1s hp:%s, hn:%s, vp:%s, vn:%s, d0:%s", r, hin, hout,
                        hout <= k ? "*" : "", toBinary(hp, m), toBinary(hn, m), toBinary(vp, m), toBinary(vn, m),
                        toBinary(d0, m));
            }

            return hout;
        }

    }

}
