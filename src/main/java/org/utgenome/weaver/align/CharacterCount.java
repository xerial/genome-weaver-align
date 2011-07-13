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
// CharacterCount.java
// Since: 2011/02/10
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;


/**
 * Character count <i>C[x]</i> is the number of characters in the input text
 * that are lexicographically smaller than the symbol <i>x</i>. This information
 * is necessary for querying texts via FM-index.
 * 
 * @author leo
 * 
 */
public class CharacterCount
{
    private static int K     = ACGT.values().length;
    private long[]     C     = new long[K];
    private long[]     count = new long[K];

    public CharacterCount(LSeq seq) {

        for (int i = 0; i < K; ++i) {
            count[i] = 0;
        }

        for (long i = 0; i < seq.textSize(); ++i) {
            ACGT ch = ACGT.decode((byte) seq.lookup(i));
            count[ch.code]++;
        }

        long sum = 0;
        for (int i = 0; i < K; ++i) {
            C[i] = sum;
            sum += count[i];
        }
    }

    public CharacterCount(WaveletArray W) {
        for (int i = 0; i < K; ++i) {
            count[i] = (int) W.rank(i, W.textSize());
        }

        long sum = 0;
        for (int i = 0; i < K; ++i) {
            C[i] = sum;
            sum += count[i];
        }
    }

    public long getCharacterCountSmallerThan(ACGT ch) {
        return C[ch.code];
    }

    public long getCount(ACGT code) {
        return count[code.code];
    }

}
