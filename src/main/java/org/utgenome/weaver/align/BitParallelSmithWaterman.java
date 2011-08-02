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

        // Searching
        for (int j = 0; j < ref.textSize(); ++j) {
            long x = pm[ref.getACGT(j).code];
            long d0 = ((vp + (x & vp)) ^ vp) | x | vn;
            long hn = vp & d0;
            long hp = (vn | ~(vp | d0));
            x = (hp << 1);
            vn = x & d0;
            vp = (hn << 1) | ~(x | d0);
            diff += (int) ((hp >>> m) & 1L);
            diff -= (int) ((hn >>> m) & 1L);
            if (_logger.isDebugEnabled()) {
                _logger.debug("j:%2d, diff:%2d %1s hp:%s, hn:%s, vp:%s, vn:%s, d0:%s", j, diff,
                        diff <= numAllowedDiff ? "*" : "", toBinary(hp, m), toBinary(hn, m), toBinary(vp, m),
                        toBinary(vn, m), toBinary(d0, m));
            }
        }
    }

    public static String toBinary(long v, int m) {
        StringBuilder s = new StringBuilder();
        for (int i = 0; i <= m; ++i) {
            s.append((v & (1L << i)) == 1L ? "1" : "0");
        }
        return s.toString();
    }

    public static void align64blocks(ACGTSequence ref, ACGTSequence query, int k) {
        final int m = (int) query.textSize();
        final int w = 64; // word size
        final int numBlocks = (int) Math.ceil((double) m / w);

        final int S = ACGT.values().length;
        // Peq bit-vector holds flags of the character occurrence positions in the query
        long[][] peq = new long[numBlocks][S];
        for (int i = 0; i < m; ++i) {
            peq[i / w][query.getACGT(i).code] |= 1L << (i % w);
        }
        // Fill the flanking region with 1s
        {
            int f = m % w;
            long mask = ~0L >>> f << f;
            for (int i = 0; i < S; ++i)
                peq[numBlocks - 1][i] |= mask;
        }

        long vp[] = new long[numBlocks];
        long vn[] = new long[numBlocks];

        final int b = (int) Math.ceil((double) k / w);
        for (int r = 0; r < b; ++r) {
            vp[r] = ~0L;
            vn[r] = 0L;
            //dt[rw * w][0] = rw; 
        }

        final int N = (int) ref.textSize();
        for (int j = 0; j < N; ++j) {
            ACGT ch = ref.getACGT(j);
            for (int r = 0; r < b; ++r) {
                int hout = alignBlock(ch, peq[r], vp[r], vn[r], w);
            }

        }

    }

    private static int alignBlock(ACGT ch, long[] peq, long vp, long vn, int w) {
        long x = peq[ch.code];
        long d0 = ((vp + (x & vp)) ^ vp) | x | vn;
        long hn = vp & d0;
        long hp = (vn | ~(vp | d0));
        x = (hp << 1);
        vn = x & d0;
        vp = (hn << 1) | ~(x | d0);

        int hout = 0;
        hout += (int) ((hp >>> w) & 1L);
        hout -= (int) ((hn >>> w) & 1L);
        return hout;
    }

    public static void align(ACGTSequence ref, ACGTSequence query) {
        final int N = (int) ref.textSize();
        final int M = (int) query.textSize();
        BitVector Pv = new BitVector(M).not();
        BitVector Mv = new BitVector(M);
        int score = M;

        int k = 2; // allowed mismatches

        // Precompute the alphabet table of query
        BitVector[] Peq = new BitVector[ACGT.values().length];
        for (int i = 0; i < Peq.length; ++i) {
            Peq[i] = new BitVector(M);
        }
        for (int m = 0; m < M; ++m) {
            Peq[query.getACGT(m).code].set(m);
        }

        for (int j = 0; j < N; ++j) {
            BitVector Eq = Peq[ref.getACGT(j).code];

            // Xv = Eq | Mv
            BitVector Xv = Eq.copy().or(Mv);
            // Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq
            BitVector Xh = Eq.copy().and(Pv).add(Pv).xor(Pv).or(Eq);

            // Ph = Mv | ~(Xv | Ph)
            // Mh = Pv & Xh
            BitVector Ph = Mv.copy().or(Xv.copy().or(Pv).not());
            BitVector Mh = Pv.copy().and(Xh);

            Ph.lshift(1);
            Mh.lshift(1);

            // 

            // Pv = Mh | ~(Xv | Ph)
            // Mv = Ph & Xv
            Pv.set(Mh.copy().or(Xv.copy().or(Ph).not()));
            Mv.set(Ph.copy().and(Xv));

            // 
            if (score <= k) {
                // match at j 
            }
        }

    }
}
