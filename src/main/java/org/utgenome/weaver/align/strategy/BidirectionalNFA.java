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
// BidirectionalNFA.java
// Since: Aug 5, 2011
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.strategy;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.PriorityQueue;

import org.utgenome.weaver.align.ACGT;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.SuffixInterval;
import org.utgenome.weaver.parallel.Reporter;
import org.xerial.util.Optional;

/**
 * NFA for forward and backward traversal of suffix arrays
 * 
 * @author leo
 * 
 */
public class BidirectionalNFA
{

    private final FMIndexOnGenome fmIndex;
    private final ACGTSequence    qF;
    private final ACGTSequence    qC;
    private final int             m;                       // query length
    private final int             numAllowedMismatches = 1;

    public BidirectionalNFA(FMIndexOnGenome fmIndex, ACGTSequence query) {
        this.fmIndex = fmIndex;
        this.qF = query;
        this.qC = query.complement();
        this.m = (int) query.textSize();
    }

    private static class CursorContainer
    {
        private List<Optional<List<Cursor>>> container;

        public CursorContainer(int m) {
            container = new ArrayList<Optional<List<Cursor>>>(m);
            for (int i = 0; i < m; ++i) {
                container.add(new Optional<List<Cursor>>());
            }
        }

        public boolean hasEntry(int index) {
            return container.get(index).isDefined();
        }

        public List<Cursor> get(int index) {
            return container.get(index).get();
        }

        public void add(int index, Cursor c) {
            Optional<List<Cursor>> l = container.get(index);
            if (l.isUndefined()) {
                l.set(new ArrayList<Cursor>());
            }
            l.get().add(c);
        }

        @Override
        public String toString() {
            return container.toString();
        }
    }

    private static enum SearchDirection {
        Forward("F"), BidirectionalForward("BF"), Backward("B");

        private final String symbol;

        private SearchDirection(String symbol) {
            this.symbol = symbol;
        }

        @Override
        public String toString() {
            return symbol;
        }
    }

    private class Cursor
    {
        public final boolean        isComplement;
        public final SuffixInterval si;
        public final int            cursorF;
        public final int            cursorB;

        public Cursor(boolean isComplement, SuffixInterval si, int cursorF, int cursorB) {
            this.isComplement = isComplement;
            this.si = si;
            this.cursorF = cursorF;
            this.cursorB = cursorB;
        }

        public int getIndex() {
            return cursorF;
        }

        public Strand getStrand() {
            return isComplement ? Strand.REVERSE : Strand.FORWARD;
        }

        public ACGT getACGT(int index) {
            return isComplement ? qC.getACGT(index) : qF.getACGT(index);
        }

        public Cursor next() {
            ACGT ch = getACGT(cursorF);
            SuffixInterval nextSi = fmIndex.forwardSearch(getStrand(), ch, si);
            if (!nextSi.isEmpty())
                return new Cursor(isComplement, nextSi, cursorF + 1, cursorB);

            // switch to bidirectional search
            return new BidirectionalCursor(isComplement, fmIndex.wholeSARange(), fmIndex.wholeSARange(), cursorF + 1,
                    cursorF + 1);
        }

        public int getRemainingBases() {
            return (m - cursorF) + cursorB;
        }

        public SearchDirection getDirection() {
            return SearchDirection.Forward;
        }

        @Override
        public String toString() {
            return String.format("%s-%d/%d:%s", getDirection(), cursorF, cursorB, si);
        }
    }

    private class BackwardCursor extends Cursor
    {
        public BackwardCursor(boolean isComplement, SuffixInterval si, int cursorF, int cursorB) {
            super(isComplement, si, cursorF, cursorB);
        }

        @Override
        public int getIndex() {
            return cursorB;
        }

        @Override
        public Cursor next() {
            ACGT ch = getACGT(cursorB - 1);
            SuffixInterval nextSi = fmIndex.backwardSearch(getStrand(), ch, si);
            if (nextSi.isEmpty())
                return null;

            return new BackwardCursor(isComplement, nextSi, cursorF, cursorB - 1);
        }

        @Override
        public SearchDirection getDirection() {
            return SearchDirection.Backward;
        }

    }

    private class BidirectionalCursor extends Cursor
    {
        public final SuffixInterval revSi;

        public BidirectionalCursor(boolean isComplement, SuffixInterval si, SuffixInterval revSi, int cursorF,
                int cursorB) {
            super(isComplement, si, cursorF, cursorB);
            this.revSi = revSi;
        }

        @Override
        public Cursor next() {
            Strand strand = getStrand();
            if (cursorF < m) {
                // Bidirectional search
                ACGT ch = getACGT(cursorF);
                SuffixInterval nextSi = fmIndex.forwardSearch(strand, ch, si);
                if (!nextSi.isEmpty()) {
                    SuffixInterval nextRevSi = fmIndex.bidirectionalSearch(strand, ch, si, revSi);
                    return new BidirectionalCursor(isComplement, nextSi, nextRevSi, cursorF + 1, cursorB);
                }

                // Start new bidirectional search
                return new BidirectionalCursor(isComplement, fmIndex.wholeSARange(), fmIndex.wholeSARange(),
                        cursorF + 1, cursorF);
            }
            else {
                // Switch to backward search
                ACGT ch = getACGT(cursorB - 1);
                SuffixInterval nextRevSi = fmIndex.backwardSearch(strand, ch, revSi);
                if (nextRevSi.isEmpty()) {
                    return null;
                }

                return new BackwardCursor(isComplement, nextRevSi, cursorF, cursorB - 1);
            }
        }

        @Override
        public SearchDirection getDirection() {
            return SearchDirection.BidirectionalForward;
        }

        @Override
        public String toString() {
            return String.format("%s-%d/%d:%s,%s", getDirection(), cursorF, cursorB, si, revSi);
        }
    }

    public static class ExactMatch
    {
        public final SuffixInterval si;

        public ExactMatch(SuffixInterval si) {
            this.si = si;
        }
    }

    /**
     * @param out
     * @throws Exception
     */
    public void bidirectionalAlign(Reporter out) throws Exception {

        // Prefer the cursor with smaller remaining bases
        PriorityQueue<Cursor> queue = new PriorityQueue<Cursor>(11, new Comparator<Cursor>() {
            @Override
            public int compare(Cursor o1, Cursor o2) {
                return o1.getRemainingBases() - o2.getRemainingBases();
            }
        });

        // Row-wise simulation of NFA
        CursorContainer stateF = new CursorContainer(m + 1);
        CursorContainer stateB = new CursorContainer(m + 1);

        Cursor initF = new Cursor(false, fmIndex.wholeSARange(), 0, 0);
        Cursor initC = new BackwardCursor(true, fmIndex.wholeSARange(), m, m);

        for (int k = 0; k < numAllowedMismatches; ++k) {
            stateF.add(initF.getIndex(), initF);
            stateB.add(initC.getIndex(), initC);
            queue.add(initF);
            queue.add(initC);

            while (!queue.isEmpty()) {
                Cursor c = queue.poll();
                if (c.getRemainingBases() == 0) {
                    // Found a match
                    continue;
                }

                Cursor next = c.next();
                if (next == null)
                    continue; // no match
                queue.add(next);

                switch (next.getDirection()) {
                case Forward:
                case BidirectionalForward:
                    stateF.add(next.getIndex(), next);
                    break;
                case Backward:
                    stateB.add(next.getIndex(), next);
                    break;
                }
            }
        }

    }

}
