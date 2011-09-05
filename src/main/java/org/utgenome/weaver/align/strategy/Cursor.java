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
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.SiSet;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.SuffixInterval;

/**
 * FM-index cursor for bidirectional search
 * 
 * @author leo
 * 
 */
public class Cursor
{
    // flag(8bit) :=  strand(1), searchDirection(2), read length(29)
    private final int flag;
    public final int  start;
    public final int  end;
    public final int  cursor;
    public final int  pivot; // switch point of forward/backward search

    private Cursor(int flag, int start, int end, int cursor, int pivot) {
        this.flag = flag;
        this.start = start;
        this.end = end;
        this.cursor = cursor;
        this.pivot = pivot;
    }

    public Cursor(Strand strand, SearchDirection searchDirection, int start, int end, int cursor, int pivot) {
        this(strand.index | (searchDirection.index << 1), start, end, cursor, pivot);
    }

    public int getFragmentLength() {
        return end - start;
    }

    @Override
    public String toString() {
        StringBuilder s = new StringBuilder();
        s.append(String.format("%s%s%d/%d[%d,%d)", getStrand().symbol, getSearchDirection().symbol, cursor, pivot,
                start, end));
        return s.toString();
    }

    public int getProcessedBases() {
        return isForwardSearch() ? cursor - pivot : end - cursor;
    }

    /**
     * @param m
     *            read length
     * @return
     */
    public int getRemainingBases() {
        return getFragmentLength() - getProcessedBases();
    }

    public Strand getStrand() {
        return Strand.decode(flag & 1);
    }

    public int getStrandIndex() {
        return flag & 1;
    }

    public int getOffsetOfSearchHead() {
        // read:   |   |------|       |
        //         0   pivot  cursor  read length
        int offset = cursor - start;
        if (getStrand() == Strand.REVERSE)
            offset = getFragmentLength() - offset;
        return offset;
    }

    public ACGT nextACGT(ACGTSequence[] q) {
        int strand = flag & 1;
        return q[strand].getACGT(getNextACGTIndex());
    }

    public int getNextACGTIndex() {
        return isForwardSearch() ? cursor : cursor - 1;
    }

    /**
     * Read search direction 0 -> m (forward), m -> 0 (backward)
     * 
     * @return
     */
    public SearchDirection getSearchDirection() {
        return SearchDirection.decode((flag >>> 1) & 0x03);
    }

    public boolean isForwardSearch() {
        return getSearchDirection().isForward;
    }

    //    public Strand fmIndexDirection() {
    //        int fm = ~(getStrandIndex() ^ (getSearchDirection().isForward ? 0 : 1)) & 1;
    //        return fm == 0 ? Strand.FORWARD : Strand.REVERSE;
    //    }

    public Cursor[] split() {
        SearchDirection d = getSearchDirection();
        Cursor left, right;
        switch (d) {
        case Forward:
            left = new Cursor(getStrand(), getSearchDirection(), start, cursor, cursor, pivot);
            right = new Cursor(getStrand(), SearchDirection.Forward, cursor, end, cursor, cursor);
            break;
        case Backward:
            left = new Cursor(getStrand(), getSearchDirection(), cursor, end, cursor, pivot);
            right = new Cursor(getStrand(), SearchDirection.Backward, start, cursor, cursor, start);
            break;
        default:
        case BidirectionalForward:
            if (cursor + 1 < end) {
                left = new Cursor(getStrand(), getSearchDirection(), start, cursor, cursor, pivot);
                right = new Cursor(getStrand(), getSearchDirection(), cursor, end, cursor, cursor);
            }
            else {
                left = new Cursor(getStrand(), SearchDirection.Backward, start, cursor, cursor, pivot);
                right = new Cursor(getStrand(), getSearchDirection(), cursor, end, cursor, cursor);
            }
            break;
        }

        return new Cursor[] { left, right };
    }

    public Cursor next() {
        int nextCursor = cursor;
        SearchDirection d = getSearchDirection();
        switch (d) {
        case Forward:
            ++nextCursor;
            break;
        case Backward:
            --nextCursor;
            break;
        case BidirectionalForward:
            if (nextCursor + 1 < end) {
                ++nextCursor;
            }
            else {
                // switch to backward search
                d = SearchDirection.Backward;
                nextCursor = pivot;
            }
            break;
        }
        return new Cursor(getStrand(), d, start, end, nextCursor, pivot);
    }

    public SiSet nextSi(FMIndexOnGenome fmIndex, SiSet si, ACGT currentBase) {
        return nextSi(fmIndex, si.getForward(currentBase), si.getBackward(currentBase));
    }

    public SiSet nextSi(FMIndexOnGenome fmIndex, SuffixInterval siF, SuffixInterval siB) {
        Strand strand = this.getStrand();
        SearchDirection d = this.getSearchDirection();
        switch (d) {
        case BidirectionalForward:
            if (cursor >= end - 1) {
                siF = null;
            }
        }
        return fmIndex.bidirectionalSearch(strand, siF, siB);
    }

}
