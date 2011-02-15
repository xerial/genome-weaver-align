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

import org.utgenome.gwt.utgb.client.bio.IUPAC;

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
    private static int K = IUPAC.values().length;
    private int[]      C = new int[K];

    public CharacterCount(LSeq seq) {

        int[] count = new int[K];
        for (int i = 0; i < K; ++i) {
            count[i] = 0;
        }

        for (long i = 0; i < seq.textSize(); ++i) {
            IUPAC code = IUPAC.decode((byte) seq.lookup(i));
            count[code.bitFlag]++;
        }
        // Decrement the character count for the sentinel since BWT based count contains one additional sentinel.   
        count[IUPAC.None.bitFlag]--;

        int sum = 1; // add 1 for $ (end of string)
        for (int i = 0; i < K; ++i) {
            C[i] = sum;
            sum += count[i];
        }

    }

    public int get(IUPAC code) {
        return C[code.bitFlag];
    }

}
