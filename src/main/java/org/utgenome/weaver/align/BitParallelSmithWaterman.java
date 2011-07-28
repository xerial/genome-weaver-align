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

/**
 * Bit-parallel algorithm for Smith-Waterman alignment
 * 
 * <h3>Reference</h3> Myers, G. (1999). A fast bit-vector algorithm for
 * approximate string matching based on dynamic programming. Journal of the ACM
 * (JACM).
 * 
 * @author leo
 * 
 */
public class BitParallelSmithWaterman
{

    public void align(ACGTSequence ref, ACGTSequence query) {
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
