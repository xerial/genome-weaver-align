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
// SiSet.java
// Since: 2011/08/12
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import org.utgenome.weaver.align.strategy.SearchDirection;

public class SiSet
{
    public final SuffixInterval[] siF;
    public final SuffixInterval[] siB;

    public SiSet(SuffixInterval[] siF, SuffixInterval[] siB) {
        this.siF = siF;
        this.siB = siB;
    }

    public SARange get(int index) {
        return new BidirectionalSuffixInterval(siF == null ? null : siF[index], siB == null ? null : siB[index]);
    }

    public boolean isEmpty(ACGT ch, SearchDirection d) {
        int i = ch.code;
        switch (d) {
        case Forward:
            return siF[i] != null && siF[i].isEmpty();
        case Backward:
            return siB[i] != null && siB[i].isEmpty();
        case BidirectionalForward:
            return (siF[i] != null && siF[i].isEmpty()) && (siB[i] != null && siB[i].isEmpty());
        }
        return false;
    }

}
