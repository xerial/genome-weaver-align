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
 * 
 * 
 * @author leo
 * 
 */
public class FMIndex
{
    private final IUPACSequence        bwt;
    private final OccurrenceCountTable O;
    private final CharacterCount       C;

    public FMIndex(IUPACSequence bwt, OccurrenceCountTable O, CharacterCount C) {
        this.bwt = bwt;
        this.O = O;
        this.C = C;
    }

    public SuffixInterval backwardSearch(IUPAC ch, SuffixInterval current) {
        int lowerBound = C.get(ch) + O.getOcc(ch, current.lowerBound - 1);
        int upperBound = C.get(ch) + O.getOcc(ch, current.upperBound) - 1;
        return new SuffixInterval(lowerBound, upperBound);
    }

    public int inverseSA(int index) {
        if (index >= bwt.size() - 1) {
            return 0;
        }
        IUPAC c = bwt.getIUPAC(index);
        return C.get(c) + O.getOcc(c, index - 1);
    }

    public int textSize() {
        return bwt.size();
    }
}
