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
// CursorBase.java
// Since: 2011/08/09
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.strategy;

import org.utgenome.weaver.align.ACGT;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.Strand;

public class CursorBase
{
    public final ACGTSequence    read;
    public final Strand          strand;
    public final SearchDirection searchDirection;
    public final int             cursorF;
    public final int             cursorB;

    public CursorBase(ACGTSequence read, Strand strand, SearchDirection searchDirection, int cursorF, int cursorB) {
        this.read = read;
        this.strand = strand;
        this.searchDirection = searchDirection;
        this.cursorF = cursorF;
        this.cursorB = cursorB;
    }

    public int getIndex() {
        return isForwardSearch() ? cursorF : cursorB - 1;
    }

    public ACGT nextACGT() {
        return isForwardSearch() ? read.getACGT(cursorF) : read.getACGT(cursorB - 1);
    }

    public boolean isForwardSearch() {
        return searchDirection.isForward;
    }

    public int getProcessedBases() {
        return cursorF - cursorB;
    }

    public int getRemainingBases() {
        return ((int) read.textSize() - cursorF) + cursorB;
    }

    public CursorBase newBidirectionalFowardCursor() {
        return new CursorBase(read, strand, SearchDirection.BidirectionalForward, cursorF + 1, cursorF + 1);
    }

    public CursorBase next() {
        SearchDirection d = searchDirection;
        int nextF = cursorF;
        int nextB = cursorB;
        switch (searchDirection) {
        case Forward:
            ++nextF;
            break;
        case Backward:
            --nextB;
            break;
        case BidirectionalForward:
            if (nextF < read.textSize()) {
                ++nextF;
            }
            else {
                // switch to backward search
                d = SearchDirection.Backward;
                --nextB;
            }
            break;
        }
        return new CursorBase(read, strand, d, nextF, nextB);
    }

}
