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
// BidirectionalCursor.java
// Since: 2011/08/08
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.strategy;

import org.utgenome.weaver.align.ACGT;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.BidirectionalSuffixInterval;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.SuffixInterval;

/**
 * cursor holding alignment state
 * 
 * @author leo
 * 
 */
public class BidirectionalCursor
{
    public final ACGTSequence    read;
    public final Strand          strand;
    public final SearchDirection searchDirection;
    public final SuffixInterval  siF;
    public final SuffixInterval  siB;
    public final int             cursorF;
    public final int             cursorB;

    public BidirectionalCursor(ACGTSequence read, Strand strand, SearchDirection searchDirection, SuffixInterval siF,
            SuffixInterval siB, int cursorF, int cursorB) {
        this.read = read;
        this.strand = strand;
        this.searchDirection = searchDirection;
        this.siF = siF;
        this.siB = siB;
        this.cursorF = cursorF;
        this.cursorB = cursorB;
    }

    public int getIndex() {
        return cursorF;
    }

    public ACGT nextACGT() {
        return isForwardSearch() ? read.getACGT(cursorF) : read.getACGT(cursorB - 1);
    }

    public boolean isForwardSearch() {
        return searchDirection.isForward;
    }

    public BidirectionalCursor next(FMIndexOnGenome fmIndex) {
        switch (searchDirection) {
        case Forward: {
            SuffixInterval nextF = fmIndex.forwardSearch(strand, nextACGT(), siF);
            if (!nextF.isEmpty())
                return new BidirectionalCursor(read, strand, searchDirection, nextF, siB, cursorF + 1, cursorB);

            // switch to bidirectional search
            return new BidirectionalCursor(read, strand, SearchDirection.BidirectionalForward, fmIndex.wholeSARange(),
                    fmIndex.wholeSARange(), cursorF + 1, cursorF + 1);
        }
        case Backward: {
            SuffixInterval nextB = fmIndex.backwardSearch(strand, nextACGT(), siB);
            if (!nextB.isEmpty())
                return new BidirectionalCursor(read, strand, searchDirection, siF, nextB, cursorF, cursorB - 1);

            // no match
            return null;
        }
        case BidirectionalForward: {
            ACGT ch = nextACGT();
            if (cursorF < (int) read.textSize()) {
                // Bidirectional search
                BidirectionalSuffixInterval next = fmIndex.bidirectionalForwardSearch(strand, ch,
                        new BidirectionalSuffixInterval(siF, siB));
                if (next != null) {
                    return new BidirectionalCursor(read, strand, searchDirection, next.forwardSi, next.backwardSi,
                            cursorF + 1, cursorB);
                }

                // Start new bidirectional search
                return new BidirectionalCursor(read, strand, searchDirection, fmIndex.wholeSARange(),
                        fmIndex.wholeSARange(), cursorF + 1, cursorF);
            }
            else {
                // Switch to backward search
                SuffixInterval nextB = fmIndex.backwardSearch(strand, ch, siB);
                if (nextB.isEmpty()) {
                    return null;
                }

                return new BidirectionalCursor(read, strand, SearchDirection.Backward, siF, nextB, cursorF, cursorB - 1);
            }
        }
        default:
            throw new IllegalStateException("cannot reach here");
        }

    }

    public int getRemainingBases() {
        return ((int) read.textSize() - cursorF) + cursorB;
    }

    @Override
    public String toString() {
        return String.format("%s%s%d/%d:%s %s", strand.symbol, searchDirection.symbol, cursorF, cursorB, siF, siB);
    }

}
