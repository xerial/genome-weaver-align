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
// FMIndex.java
// Since: 2011/02/12
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import org.utgenome.gwt.utgb.client.bio.IUPAC;

/**
 * Text-search index based on BWT string
 * 
 * @author leo
 * 
 */
public interface FMIndex
{
    public SuffixInterval backwardSearch(IUPAC ch, SuffixInterval current);

    public CharacterCount getCharacterCount();

    /**
     * Follow the suffix link using the equation: SA[x] - 1 = C(x) + Rank(c, x).
     * 
     * @param index
     *            index x in the suffix array
     * @return index p in the suffix array that satisfies SA[p] = SA[x] - 1.
     */
    public long suffixLink(long index);

    public long textSize();
}
