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
import org.utgenome.weaver.align.AlignmentScoreConfig;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.SuffixInterval;

/**
 * cursor holding an alignment state of the bidirectional search
 * 
 * @author leo
 * 
 */
public class BidirectionalCursor
{
    public final Score               score;
    public final ACGTSequence        read;
    public final Strand              strand;
    public final SearchDirection     searchDirection;
    public final ExtensionType       extensionType;
    public final SuffixInterval      siF;
    public final SuffixInterval      siB;
    public final int                 cursorF;
    public final int                 cursorB;
    public final BidirectionalCursor split;

    public BidirectionalCursor(Score score, ACGTSequence read, Strand strand, SearchDirection searchDirection,
            ExtensionType extensionType, SuffixInterval siF, SuffixInterval siB, int cursorF, int cursorB,
            BidirectionalCursor split) {
        this.score = score;
        this.read = read;
        this.strand = strand;
        this.searchDirection = searchDirection;
        this.extensionType = extensionType;
        this.siF = siF;
        this.siB = siB;
        this.cursorF = cursorF;
        this.cursorB = cursorB;
        this.split = split;
    }

    public BidirectionalCursor(Score score, ACGTSequence read, Strand strand, SearchDirection searchDirection,
            ExtensionType extensionType, SuffixInterval siF, SuffixInterval siB, int cursorF, int cursorB) {
        this(score, read, strand, searchDirection, extensionType, siF, siB, cursorF, cursorB, null);
    }

    public int getUpperBoundOfScore(AlignmentScoreConfig config) {
        return score.score + getRemainingBases() * config.matchScore;
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

    public int getRemainingBases() {
        return ((int) read.textSize() - cursorF) + cursorB;
    }

    @Override
    public String toString() {
        StringBuilder s = new StringBuilder();
        s.append(String.format("%s%s%2d:%d/%d:%s", strand.symbol, searchDirection.symbol, score.score, cursorF,
                cursorB, siF));
        if (siB != null)
            s.append(String.format(" %s", siB));
        return s.toString();
    }

    //public BidirectionalCursor extendWithMatch(AlignmentScoreConfig config, )

}
