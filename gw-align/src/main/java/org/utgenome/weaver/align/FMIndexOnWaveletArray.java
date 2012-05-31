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
// FMIndexOnWaveletArray.java
// Since: 2011/05/16
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

public class FMIndexOnWaveletArray implements FMIndex
{
    private final WaveletArray   W;
    private final CharacterCount C;

    public FMIndexOnWaveletArray(WaveletArray W) {
        this.W = W;
        this.C = new CharacterCount(W);
    }

    public CharacterCount getCharacterCount() {
        return C;
    }

    public SuffixInterval backwardSearch(ACGT ch, SuffixInterval current) {
        long lowerBound = C.getCharacterCountSmallerThan(ch) + W.rank(ch.code, current.lowerBound);
        long upperBound = C.getCharacterCountSmallerThan(ch) + W.rank(ch.code, current.upperBound);
        return new SuffixInterval(lowerBound, upperBound);
    }

    /**
     * Follow the suffix link using the equation: SA[x] - 1 = C(x) + Rank(c, x).
     * 
     * @param index
     *            index x in the suffix array
     * @return index p in the suffix array that satisfies SA[p] = SA[x] - 1.
     */
    public long suffixLink(long index) {
        if (index >= W.textSize()) { // If the index reaches the sentinel 
            return 0; // Return the smallest SA index
        }
        ACGT c = ACGT.decode((byte) W.lookup(index));
        return C.getCharacterCountSmallerThan(c) + W.rank(c.code, index);
    }

    public long textSize() {
        return W.textSize();
    }

    @Override
    public long count(ACGT ch, long start, long end) {
        throw new UnsupportedOperationException();
    }

    @Override
    public long[] rankACGTN(long suffixIndex) {
        long[] rank = new long[ACGT.values().length];
        for (ACGT ch : ACGT.values())
            rank[ch.code] = W.rank(ch.code, suffixIndex);
        return rank;
    }
}
