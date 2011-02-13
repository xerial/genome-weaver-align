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
// SparseSuffixArray.java
// Since: 2011/02/12
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

/**
 * 
 * 
 * S(i) = S(X^{(j)}(i)) + j
 * 
 * X(i) = C(BWT[i]) + O(BWT[i], i)
 * 
 * @author leo
 * 
 */
public class SparseSuffixArray
{
    private final int[]   sparseSA;
    private final int     L;
    private final FMIndex fmIndex;

    public SparseSuffixArray(int[] sparseSA, int L, FMIndex fmIndex) {
        this.sparseSA = sparseSA;
        this.fmIndex = fmIndex;
        this.L = L;
    }

    public int get(int index) {
        int pos = index / L;
        int offset = index % L;
        if (offset == 0)
            return sparseSA[pos];

        int cursor = index;
        final int N = fmIndex.size();
        for (int j = 0; j < N; j++) {
            cursor = fmIndex.inverseSA(cursor);
            if (cursor % L == 0)
                return sparseSA[cursor / L] + j;
        }
        throw new IllegalStateException("cannot reach here");

    }
}
