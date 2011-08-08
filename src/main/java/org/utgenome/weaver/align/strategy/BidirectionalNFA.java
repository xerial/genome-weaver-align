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
import org.utgenome.weaver.align.BidirectionalSuffixInterval;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.SuffixInterval;
import org.utgenome.weaver.parallel.Reporter;
import org.xerial.lens.SilkLens;
import org.xerial.util.Optional;
import org.xerial.util.log.Logger;

/**
 * NFA for forward and backward traversal of suffix arrays
 * 
 * @author leo
 * 
 */
public class BidirectionalNFA
{
    private static Logger         _logger              = Logger.getLogger(BidirectionalNFA.class);

    private final FMIndexOnGenome fmIndex;
    private final ACGTSequence    qF;
    private final ACGTSequence    qC;
    private final int             m;                                                              // query length
    private final int             numAllowedMismatches = 1;

    public BidirectionalNFA(FMIndexOnGenome fmIndex, ACGTSequence query) {
        this.fmIndex = fmIndex;
        this.qF = query;
        this.qC = query.complement();
        this.m = (int) query.textSize();
    }

    private static class CursorContainer
    {
        private List<Optional<List<Cursor>>> containerF;
        private List<Optional<List<Cursor>>> containerR;

        public CursorContainer(int m) {
            containerF = new ArrayList<Optional<List<Cursor>>>(m);
            containerR = new ArrayList<Optional<List<Cursor>>>(m);
            for (int i = 0; i < m; ++i) {
                containerF.add(new Optional<List<Cursor>>());
                containerR.add(new Optional<List<Cursor>>());
            }
        }

        private List<Optional<List<Cursor>>> getContainer(Strand strand) {
            return strand.isForward() ? containerF : containerR;
        }

        public void add(Cursor c) {
            Optional<List<Cursor>> l = getContainer(c.strand).get(c.getIndex());
            if (l.isUndefined()) {
                l.set(new ArrayList<Cursor>());
            }
            l.get().add(c);
        }

        @Override
        public String toString() {
            return String.format("%s\n%s", containerF, containerR);
        }
    }

    private abstract class Cursor
    {
        public final Strand         strand;
        public final SuffixInterval si;
        public final int            cursorF;
        public final int            cursorB;

        public Cursor(Strand strand, SuffixInterval si, int cursorF, int cursorB) {
            this.strand = strand;
            this.si = si;
            this.cursorF = cursorF;
            this.cursorB = cursorB;
        }

        public int getIndex() {
            return cursorF;
        }

        public ACGT nextACGT() {
            ACGTSequence q = strand.isForward() ? qF : qC;
            return isForwardSearch() ? q.getACGT(cursorF) : q.getACGT(cursorB - 1);
        }

        public boolean isForwardSearch() {
            return true;
        }

        public abstract String getSymbol();

        public abstract Cursor next();

        public int getRemainingBases() {
            return (m - cursorF) + cursorB;
        }

        @Override
        public String toString() {
            return String.format("%s%s%d/%d:%s", this.getSymbol(), strand.symbol, cursorF, cursorB, si);
        }
    }

    private class ForwardCursor extends Cursor
    {
        public ForwardCursor(Strand strand, SuffixInterval si, int cursorF, int cursorB) {
            super(strand, si, cursorF, cursorB);
        }

        @Override
        public Cursor next() {
            SuffixInterval nextSi = fmIndex.forwardSearch(strand, nextACGT(), si);
            if (!nextSi.isEmpty())
                return new ForwardCursor(strand, nextSi, cursorF + 1, cursorB);

            // switch to bidirectional search
            return new BidirectionalForwardCursor(strand, fmIndex.wholeSARange(), fmIndex.wholeSARange(), cursorF + 1,
                    cursorF + 1);
        }

        @Override
        public String getSymbol() {
            return "F";
        }

    }

    private class BackwardCursor extends Cursor
    {
        public BackwardCursor(Strand strand, SuffixInterval si, int cursorF, int cursorB) {
            super(strand, si, cursorF, cursorB);
        }

        @Override
        public int getIndex() {
            return cursorB;
        }

        @Override
        public Cursor next() {
            SuffixInterval nextSi = fmIndex.backwardSearch(strand, nextACGT(), si);
            if (!nextSi.isEmpty())
                return new BackwardCursor(strand, nextSi, cursorF, cursorB - 1);

            // no match
            return null;
        }

        @Override
        public boolean isForwardSearch() {
            return false;
        }

        @Override
        public String getSymbol() {
            return "R";
        }
    }

    private class BidirectionalForwardCursor extends Cursor
    {
        public final SuffixInterval revSi;

        public BidirectionalForwardCursor(Strand readOrientation, SuffixInterval si, SuffixInterval revSi, int cursorF,
                int cursorB) {
            super(readOrientation, si, cursorF, cursorB);
            this.revSi = revSi;
        }

        @Override
        public Cursor next() {
            ACGT ch = nextACGT();
            if (cursorF < m) {
                // Bidirectional search
                BidirectionalSuffixInterval next = fmIndex.bidirectionalForwardSearch(strand, ch,
                        new BidirectionalSuffixInterval(si, revSi));
                if (next != null) {
                    return new BidirectionalForwardCursor(strand, next.forwardSi, next.backwardSi, cursorF + 1, cursorB);
                }

                // Start new bidirectional search
                return new BidirectionalForwardCursor(strand, fmIndex.wholeSARange(), fmIndex.wholeSARange(),
                        cursorF + 1, cursorF);
            }
            else {
                // Switch to backward search
                SuffixInterval nextRevSi = fmIndex.backwardSearch(strand, ch, revSi);
                if (nextRevSi.isEmpty()) {
                    return null;
                }

                return new BackwardCursor(strand, nextRevSi, cursorF, cursorB - 1);
            }
        }

        @Override
        public String toString() {
            return String.format("%s%s%d/%d:%s,%s", getSymbol(), strand.symbol, cursorF, cursorB, si, revSi);
        }

        @Override
        public String getSymbol() {
            return "B";
        }
    }

    public static class ExactMatch
    {
        public final Strand         strand;
        public final SuffixInterval si;
        public final int            numMismatches;

        public ExactMatch(Strand strand, SuffixInterval si, int numMismatches) {
            this.strand = strand;
            this.si = si;
            this.numMismatches = numMismatches;
        }
    }

    public void align() throws Exception {
        bidirectionalAlign(new Reporter() {
            @Override
            public void emit(Object result) throws Exception {
                _logger.debug(SilkLens.toSilk(result));
            }
        });
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
        CursorContainer state = new CursorContainer(m + 1);

        Cursor initF = new ForwardCursor(Strand.FORWARD, fmIndex.wholeSARange(), 0, 0);
        Cursor initC = new ForwardCursor(Strand.REVERSE, fmIndex.wholeSARange(), 0, 0);
        queue.add(initF);
        queue.add(initC);
        state.add(initF);
        state.add(initC);

        for (int k = 0; k < numAllowedMismatches; ++k) {

            while (!queue.isEmpty()) {
                Cursor c = queue.poll();
                if (c.getRemainingBases() == 0) {
                    // Found a match
                    out.emit(new ExactMatch(c.strand, c.si, k));
                    continue;
                }

                Cursor next = c.next();
                if (next == null)
                    continue; // no match
                queue.add(next);

                state.add(next);
            }
        }

    }

}
