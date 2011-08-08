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
    private final ACGTSequence    query;
    private final Strand          strand;
    private final int             m;      // query length

    public BidirectionalNFA(FMIndexOnGenome fmIndex, ACGTSequence query, Strand strand) {
        this.fmIndex = fmIndex;
        this.query = query;
        this.strand = strand;
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
        public final SuffixInterval si;
        public final int            cursorF;
        public final int            cursorB;

        public Cursor(SuffixInterval si, int cursorF, int cursorB) {
            this.si = si;
            this.cursorF = cursorF;
            this.cursorB = cursorB;
        }

        public int getIndex() {
            return cursorF;
        }

        public Cursor next() {
            ACGT ch = query.getACGT(cursorF);
            SuffixInterval nextSi = fmIndex.forwardSearch(strand, ch, si);
            if (!nextSi.isEmpty())
                return new Cursor(nextSi, cursorF + 1, cursorB);

            // switch to bidirectional search
            return new BidirectionalCursor(fmIndex.wholeSARange(), fmIndex.wholeSARange(), cursorF + 1, cursorF + 1);
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
        public BackwardCursor(SuffixInterval si, int cursorF, int cursorB) {
            super(si, cursorF, cursorB);
        }

        @Override
        public int getIndex() {
            return cursorB;
        }

        @Override
        public Cursor next() {
            ACGT ch = query.getACGT(cursorB);
            SuffixInterval nextSi = fmIndex.backwardSearch(ch, si);
            if (nextSi.isEmpty())
                return null;

            return new BackwardCursor(nextSi, cursorF, cursorB - 1);
        }

        @Override
        public SearchDirection getDirection() {
            return SearchDirection.Backward;
        }

    }

    private class BidirectionalCursor extends Cursor
    {
        public final SuffixInterval revSi;

        public BidirectionalCursor(SuffixInterval si, SuffixInterval revSi, int cursorF, int cursorB) {
            super(si, cursorF, cursorB);
            this.revSi = revSi;
        }

        @Override
        public Cursor next() {
            if (cursorF < m) {
                // Bidirectional search
                ACGT ch = query.getACGT(cursorF);
                SuffixInterval nextSi = fmIndex.forwardSearch(strand, ch, si);
                if (!nextSi.isEmpty()) {
                    SuffixInterval nextRevSi = fmIndex.bidirectionalSearch(strand, ch, si, revSi);
                    return new BidirectionalCursor(nextSi, nextRevSi, cursorF + 1, cursorB);
                }

                // Start new bidirectional search
                return new BidirectionalCursor(fmIndex.wholeSARange(), fmIndex.wholeSARange(), cursorF + 1, cursorF);
            }
            else {
                // Switch to backward search
                ACGT ch = query.getACGT(cursorB);
                SuffixInterval nextRevSi = fmIndex.backwardSearch(ch, revSi);
                if (nextRevSi.isEmpty()) {
                    return null;
                }

                return new BackwardCursor(nextRevSi, cursorF, cursorB - 1);
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

        Cursor init = new Cursor(fmIndex.wholeSARange(), 0, 0);
        stateF.add(init.getIndex(), init);
        queue.add(init);

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

    public void align(Reporter out) throws Exception {

        // Row-wise simulation of NFA
        CursorContainer stateF = new CursorContainer(m + 1);
        CursorContainer stateB = new CursorContainer(m + 1);

        // Simulate k=0 (no mismatch)
        // Forward search
        {
            SuffixInterval si = fmIndex.wholeSARange();
            int breakPoint = 0;
            stateF.add(0, new Cursor(si, 0, 0));
            for (int i = 0; i < m; ++i) {
                ACGT ch = query.getACGT(i);
                SuffixInterval nextSi = fmIndex.forwardSearch(strand, ch, si);
                if (!nextSi.isEmpty()) {
                    si = nextSi;
                }
                else {
                    breakPoint = i;
                    si = fmIndex.wholeSARange();
                }
                stateF.add(i + 1, new Cursor(si, i, breakPoint));
            }

            if (breakPoint == 0) {
                // exact match
                out.emit(new ExactMatch(si));
            }

        }

        // Backward search from the last match
        if (stateF.hasEntry(m)) {
            for (Cursor each : stateF.get(m)) {
                // compute reverse SA range
                SuffixInterval revSi = each.si;
                for (int i = each.cursorB; i < each.cursorF; ++i) {
                    ACGT ch = query.getACGT(i);
                    revSi = fmIndex.bidirectionalSearch(strand, ch, each.si, revSi);
                }
                stateB.add(each.cursorB, new Cursor(revSi, each.cursorB, each.cursorF));
            }
        }

    }

}
