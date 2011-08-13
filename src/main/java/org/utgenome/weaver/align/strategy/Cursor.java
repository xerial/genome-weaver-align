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
// Cursor.java
// Since: Aug 13, 2011
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.strategy;

import org.utgenome.weaver.align.ACGT;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.Strand;

/**
 * FM-index cursor for bidirectional search
 * 
 * @author leo
 * 
 */
public class Cursor
{
    // flag(8bit) :=  strand(1), searchDirection(2), read length(29)
    private final int  flag;
    public final int    cursorF;
    public final int    cursorB;
    public final Cursor split;

    private Cursor(int flag, int cursorF, int cursorB, Cursor split) {
        this.flag = flag;
        this.cursorF = cursorF;
        this.cursorB = cursorB;
        this.split = split;
    }

    public Cursor(Strand strand, SearchDirection searchDirection, int readLength, int cursorF, int cursorB, Cursor split) {
        this(strand.index | (searchDirection.index << 1) | (readLength << 3), cursorF, cursorB, split);
    }

    public int getReadLength() {
    	return flag >>> 3;
    }
    
    @Override
    public String toString() {
        StringBuilder s = new StringBuilder();
        s.append(String.format("%s%s:%d/%d", getStrand().symbol, getSearchDirection().symbol, cursorF, cursorB));
        if (split != null)
            s.append(String.format(" split(%s)", split));
        return s.toString();
    }

    public int getProcessedBases() {
        int p = getCursorRange();
        if (split != null)
            p += split.getCursorRange();
        return p;
    }

    private int getCursorRange() {
    	return cursorF - cursorB;
    }
    
    /**
     * @param m	read length
     * @return
     */
    public int getRemainingBases() {
    	return getReadLength() - getProcessedBases();
    }

    public Strand getStrand() {
        return Strand.decode(flag & 1);
    }

    public int getStrandIndex() {
        return flag & 1;
    }

    public ACGT nextACGT(ACGTSequence[] q) {
        int strand = flag & 1;
        return q[strand].getACGT(getNextACGTIndex());
    }

    public int getNextACGTIndex() {
        return getSearchDirection().isForward && cursorF < getReadLength() ? cursorF : cursorB - 1;
    }

    public int getIndex() {
        if (split == null)
            return cursorF - cursorB;
        else
            return cursorF - cursorB + (split.cursorF - split.cursorB);
    }

    public SearchDirection getSearchDirection() {
        return SearchDirection.decode((flag >>> 1) & 0x03);
    }

    Cursor split() {
        int cursor = cursorF;
        if (getSearchDirection() == SearchDirection.Backward)
            cursor = cursorB;

        return new Cursor(flag, cursor, cursor, this);
    }

    public boolean hasSplit() {
        return split != null;
    }

    public Cursor next() {
        Strand strand = getStrand();
        int nextF = cursorF;
        int nextB = cursorB;
        SearchDirection d = getSearchDirection();
        switch (d) {
        case Forward:
            ++nextF;
            break;
        case Backward:
            --nextB;
            break;
        case BidirectionalForward:
            if (nextF < getReadLength()) {
                ++nextF;
            }
            else {
                // switch to backward search
                d = SearchDirection.Backward;
                --nextB;
            }
            break;
        }
        return new Cursor(strand, d, getReadLength(), nextF, nextB, split);
    }
}



