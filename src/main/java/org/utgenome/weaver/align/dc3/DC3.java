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
// DC3.java
// Since: 2011/04/07
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.dc3;

import java.util.ArrayList;
import java.util.Random;

import org.utgenome.weaver.align.LSeq;
import org.utgenome.weaver.align.sais.UInt32Array;
import org.xerial.util.StringUtil;

/**
 * DC3 Suffix array building algorithm
 * 
 * @author leo
 * 
 */
public class DC3
{

    public static UInt32Array buildSuffixArray(LSeq T) {

        final long N = (T.textSize() * 2 / 3) + (T.textSize() % 3);

        // S := {(T[i, i+2], i): i \in [0, N), i mod 3 != 0 }
        UInt32Array S = new UInt32Array(N);
        for (long i = 0, j = 0; i < T.textSize(); ++i) {
            if (i % 3 == 0)
                continue;
            S.set(j++, i);
        }

        // sort S by the 3-mer T[S[i], S[i]+2]
        sortTrimer(T, S);

        return null;
    }

    public static void sortTrimer(LSeq T, UInt32Array S) {
        TrimerQuickSort.sort(T, S);
    }

    public static class TrimerQuickSort
    {
        private final LSeq        T;
        private final UInt32Array S;
        private final Random      random = new Random(0);

        @Override
        public String toString() {
            ArrayList<String> l = new ArrayList<String>();
            for (long i = 0; i < S.textSize(); ++i) {
                long idx = S.lookup(i);
                StringBuilder seq = new StringBuilder(3);
                for (long x = idx; x < idx + 3 && x < T.textSize(); ++x) {
                    seq.append((char) T.lookup(x));
                }
                l.add(String.format("(%s, %d)", seq.toString(), idx));
            }
            return StringUtil.join(l, ", ");
        }

        private TrimerQuickSort(LSeq T, UInt32Array S) {
            this.T = T;
            this.S = S;
        }

        public static void sort(LSeq T, UInt32Array S) {
            new TrimerQuickSort(T, S).quickSort(0, S.textSize() - 1);
        }

        private void swap(long i, long j) {
            long tmp = S.lookup(i);
            S.set(i, S.lookup(j));
            S.set(j, tmp);
        }

        private void quickSort(long lowerBound, long upperBound) {
            if (lowerBound >= upperBound)
                return;

            long pivotIndex = (Math.abs(random.nextLong()) % (upperBound - lowerBound)) + lowerBound;
            long pivot = S.lookup(pivotIndex);
            // move pivot to the last of the array 
            swap(pivotIndex, upperBound);

            // partition
            long splitIndex = lowerBound;
            for (long i = lowerBound; i < upperBound - 1; ++i) {
                long x = S.lookup(i);
                if (compareTrimer(x, pivot) < 0) {
                    swap(splitIndex, i);
                    splitIndex++;
                }
            }
            swap(splitIndex, upperBound);

            quickSort(lowerBound, splitIndex - 1);
            quickSort(splitIndex + 1, upperBound);
        }

        private int compareTrimer(long x, long y) {
            long N = T.textSize();
            for (int i = 0; i < 3; ++i) {
                long Sx = x + i < N ? T.lookup(x + i) : 0;
                long Sy = y + i < N ? T.lookup(y + i) : 0;
                if (Sx == Sy)
                    continue;
                else
                    return (Sx < Sy) ? -1 : 1;
            }
            return 0;
        }

    }

}
