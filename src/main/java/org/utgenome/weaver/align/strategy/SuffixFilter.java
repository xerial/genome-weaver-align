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
// SuffixFilter.java
// Since: 2011/07/27
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.strategy;

import java.util.Comparator;
import java.util.PriorityQueue;

import org.utgenome.weaver.align.ACGT;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.AlignmentScoreConfig;
import org.utgenome.weaver.align.BidirectionalSuffixInterval;
import org.utgenome.weaver.align.BitVector;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.SuffixInterval;
import org.utgenome.weaver.parallel.Reporter;
import org.xerial.lens.SilkLens;
import org.xerial.util.log.Logger;

/**
 * Suffix filter
 * 
 * <pre>
 *   *---*---*---*
 *   | \ | \ | \ |
 *   *---*---*---*
 *   | \ | \ | \ |
 *   *---*---*---*
 * 
 * </pre>
 * 
 * @author leo
 * 
 */
public class SuffixFilter
{
    private static Logger              _logger = Logger.getLogger(SuffixFilter.class);

    private final FMIndexOnGenome      fmIndex;
    private final AlignmentScoreConfig config;
    private final int                  k;
    private final int                  m;

    private StaircaseFilter            staircaseFilter;

    public static enum Diff {
        Match, Mismatch, Insertion, Deletion, Split
    };

    /**
     * FM-index cursor for bidirectional search
     * 
     * @author leo
     * 
     */
    public class Cursor
    {
        // flag(8bit) :=  strand(1), searchDirection(2)
        private final byte          flag;
        //private Score      score;
        public final SuffixInterval siF;
        public final SuffixInterval siB;
        public final int            cursorF;
        public final int            cursorB;
        public final Cursor         split;

        private Cursor(byte flag, SuffixInterval siF, SuffixInterval siB, int cursorF, int cursorB, Cursor split) {
            this.flag = flag;
            this.siF = siF;
            this.siB = siB;
            this.cursorF = cursorF;
            this.cursorB = cursorB;
            this.split = split;
        }

        public Cursor(Strand strand, SearchDirection searchDirection, SuffixInterval siF, SuffixInterval siB,
                int cursorF, int cursorB, Cursor split) {
            this((byte) (strand.index | (searchDirection.index << 1)), siF, siB, cursorF, cursorB, split);
        }

        @Override
        public String toString() {
            StringBuilder s = new StringBuilder();
            s.append(String.format("%s%s:%d/%d:%s", getStrand().symbol, getSearchDirection().symbol, cursorF, cursorB,
                    siF));
            if (siB != null)
                s.append(String.format(" %s", siB));
            if (split != null)
                s.append(String.format(" split(%s)", split));
            return s.toString();
        }

        public int getProcessedBases() {
            int p = cursorF - cursorB;
            if (split == null)
                return p;
            else
                return p + (split.cursorF - split.cursorB);
        }

        public int getRemainingBases() {
            int r = (m - cursorF) + cursorB;
            if (split == null)
                return r;
            else
                return r - (split.cursorF - split.cursorB);
        }

        public Strand getStrand() {
            return Strand.decode(flag & 1);
        }

        public ACGT nextACGT(ACGTSequence[] q) {
            int strand = flag & 1;
            return getSearchDirection().isForward && cursorF < m ? q[strand].getACGT(cursorF) : q[strand]
                    .getACGT(cursorB - 1);
        }

        public int getIndex() {
            return cursorF - cursorB;
        }

        public SearchDirection getSearchDirection() {
            return SearchDirection.decode((flag >>> 1) & 0x03);
        }

        public ExtensionType getExtentionType() {
            return ExtensionType.decode((flag >>> 3) & 0x03);
        }

        public Cursor split() {
            int cursor = cursorF;
            if (getSearchDirection() == SearchDirection.Backward)
                cursor = cursorB;

            return new Cursor(flag, fmIndex.wholeSARange(), fmIndex.wholeSARange(), cursor, cursor, this);
        }

        public Cursor next(ACGT ch) {
            Strand strand = getStrand();
            SuffixInterval nextSiF = siF;
            SuffixInterval nextSiB = siB;
            int nextF = cursorF;
            int nextB = cursorB;
            SearchDirection d = getSearchDirection();
            switch (d) {
            case Forward:
                nextSiF = fmIndex.forwardSearch(strand, ch, siF);
                if (nextSiF.isEmpty())
                    return null;
                ++nextF;
                break;
            case Backward:
                nextSiB = fmIndex.backwardSearch(strand, ch, siB);
                if (nextSiB.isEmpty())
                    return null;
                --nextB;
                break;
            case BidirectionalForward:
                if (nextF < m) {
                    ++nextF;
                    BidirectionalSuffixInterval bSi = fmIndex.bidirectionalForwardSearch(strand, ch,
                            new BidirectionalSuffixInterval(siF, siB));
                    if (bSi == null)
                        return null;
                    nextSiF = bSi.forwardSi;
                    nextSiB = bSi.backwardSi;
                }
                else {
                    // switch to backward search
                    d = SearchDirection.Backward;
                    --nextB;
                    nextSiB = fmIndex.backwardSearch(strand, ch, siB);
                    if (nextSiB.isEmpty())
                        return null;
                }
                break;
            }

            return new Cursor(flag, nextSiF, nextSiB, nextF, nextB, split);
        }
    }

    /**
     * NFA state and alignment cursor holder
     * 
     * TODO optimize the automaton size
     * 
     * @author leo
     * 
     */
    public class SearchState
    {
        public final Cursor       cursor;
        private final BitVector[] automaton;
        // 32 bit = searchFlag (5) + minK (8)
        private int               state;

        private SearchState(Cursor cursor, BitVector[] automaton, int minK) {
            this.cursor = cursor;
            this.automaton = automaton;
            this.state = minK << 5;
        }

        @Override
        public String toString() {
            return String.format("k%d%s", getNumDifferences(), cursor);
        }

        public int getStrandIndex() {
            return cursor.flag & 1;
        }

        public int getIndex() {
            return cursor.getIndex();
        }

        public int getNumDifferences() {
            return (state >>> 5) & 0xFF;
        }

        public boolean isFinished() {
            return (state & 0x1F) == 0x1F; // ACGT + split
        }

        public void updateFlag(ACGT ch) {
            this.state |= 1 << ch.code;
        }

        public void updateSplitFlag() {
            this.state |= 1 << 4;
        }

        public boolean isChecked(ACGT ch) {
            return (state & (1 << ch.code)) != 0;
        }

        public SearchState nextStateAfterSplit() {
            updateSplitFlag();
            // use the same automaton state
            if (getNumDifferences() < k) {
                return new SearchState(cursor.split(), automaton, getNumDifferences() + 1);
            }
            else
                return null;
        }

        public SearchState nextState(ACGT ch, Cursor nextCursor, QueryMask queryMask) {

            BitVector[] prev = automaton;
            BitVector[] next = new BitVector[k + 1];

            final BitVector qeq = queryMask.getPatternMask(ch);

            final int nextIndex = cursor.getIndex() + 1;

            int nm = 0;
            // R'_0 = (R_0 << 1) | P[ch]
            next[nm] = prev[nm].rshift(1)._and(qeq);
            if (next[nm].get(m - 1)) {
                // Found a full match
                //return new DFSState(State.Finished, 0);
                return new SearchState(nextCursor, next, nm);
            }
            if (!next[nm].get(nextIndex))
                ++nm;

            for (int i = 1; i <= k; ++i) {
                // R'_{i+1} = ((R_{i+1} << 1) &  P[ch]) | R_i | (R_i << 1) | (R'_i << 1)   
                next[i] = prev[i].rshift(1)._and(qeq);
                next[i]._or(prev[i - 1]);
                next[i]._or(prev[i - 1].rshift(1));
                next[i]._or(next[i - 1].rshift(1));
                // Apply a suffix filter (staircase mask)
                next[i]._and(staircaseFilter.getStaircaseMask(i));

                // Found a match
                if (next[i].get(m - 1))
                    return new SearchState(nextCursor, next, i);

                if (nm == i && !next[i].get(nextIndex))
                    ++nm;
            }

            if (nm >= k) {
                // no match
                return null;
            }
            else {
                // extend the match
                return new SearchState(nextCursor, next, nm);
            }
        }

    }

    private SearchState initialState(Strand strand, SearchDirection searchDirection, int k, int m) {
        BitVector[] automaton = new BitVector[k + 1];
        // Activate the diagonal states 
        for (int i = 0; i <= k; ++i) {
            automaton[i] = new BitVector(m);
            // TODO optimize flag set
            for (int j = 0; j <= i; ++j)
                automaton[i].set(j);
        }
        SearchState s = new SearchState(new Cursor(strand, searchDirection, fmIndex.wholeSARange(), null, 0, 0, null),
                automaton, 0);

        return s;
    }

    /**
     * A set of bit flags of ACGT characters in a query sequence
     * 
     * @author leo
     * 
     */
    public static class QueryMask
    {
        private BitVector[] patternMask;

        public QueryMask(ACGTSequence query) {
            this(query, 0);
        }

        public QueryMask(ACGTSequence query, int offset) {
            int m = (int) query.textSize();
            patternMask = new BitVector[ACGT.exceptN.length];
            for (int i = 0; i < patternMask.length; ++i)
                patternMask[i] = new BitVector(m + 1);

            for (int i = 0; i < m; ++i) {
                int index = offset + i;
                if (index >= m) {
                    // for bidirectional search
                    index = m - i - 1;
                }
                ACGT ch = query.getACGT(index);
                if (ch == ACGT.N) {
                    for (ACGT each : ACGT.exceptN)
                        patternMask[each.code].set(index + 1);
                }
                else
                    patternMask[ch.code].set(index + 1);
            }
        }

        public BitVector getPatternMask(ACGT ch) {
            return patternMask[ch.code];
        }
    }

    public SuffixFilter(FMIndexOnGenome fmIndex, AlignmentScoreConfig config, long m) {
        this.fmIndex = fmIndex;
        this.config = config;
        this.k = config.maximumEditDistances;
        this.m = (int) m;
        this.staircaseFilter = new StaircaseFilter(this.m, k);
    }

    public void align(ACGTSequence query) throws Exception {
        new AlignmentProcess(query, new Reporter() {
            @Override
            public void emit(Object result) throws Exception {
                _logger.debug(SilkLens.toSilk("result", result));

            }
        }).align();
    }

    public void align(ACGTSequence query, Reporter out) throws Exception {
        new AlignmentProcess(query, out).align();
    }

    private static class StatePreference implements Comparator<SearchState>
    {
        @Override
        public int compare(SearchState o1, SearchState o2) {
            // prefer longer match
            int diff = o2.getIndex() - o1.getIndex();
            if (diff != 0)
                return diff;

            int kDiff = o1.getNumDifferences() - o2.getNumDifferences();
            return kDiff;
        }
    }

    public class AlignmentProcess
    {

        private ACGTSequence[]             q         = new ACGTSequence[2];
        private QueryMask[]                queryMask = new QueryMask[2];
        private PriorityQueue<SearchState> queue     = new PriorityQueue<SearchState>(11, new StatePreference());
        private Reporter                   out;

        public AlignmentProcess(ACGTSequence query, Reporter out) {
            this.q[0] = query;
            this.q[1] = query.complement();
            this.queryMask[0] = new QueryMask(q[0]);
            this.queryMask[1] = new QueryMask(q[1]);
            this.out = out;
        }

        public void align() throws Exception {

            if (q[0].fastCount(ACGT.N, 0, m) > config.maximumEditDistances) {
                return; // skip this alignment
            }

            // Add states for both strands
            queue.add(initialState(Strand.FORWARD, SearchDirection.Forward, k, m));
            queue.add(initialState(Strand.REVERSE, SearchDirection.Forward, k, m));

            while (!queue.isEmpty()) {
                SearchState c = queue.poll();
                if (c.isFinished()) {
                    continue;
                }

                int strandIndex = c.getStrandIndex();
                int allowedMismatches = k - c.getNumDifferences();
                if (c.getIndex() + allowedMismatches >= m) {
                    // found a hit
                    reportAlignment(c);
                    continue;
                }

                if (allowedMismatches == 0) {
                    // do exact match
                    SearchState matchState = exactMatch(c);
                    if (matchState != null)
                        reportAlignment(matchState);
                    continue;
                }

                ACGT nextBase = c.cursor.nextACGT(q);
                if (!c.isChecked(nextBase)) {
                    // search for a base in the read
                    step(c, nextBase);
                    queue.add(c); // preserve the state for backtracking
                    continue;
                }

                // search for the other bases
                for (ACGT ch : ACGT.exceptN) {
                    if (ch == nextBase)
                        continue;

                    if (!c.isChecked(ch)) {
                        step(c, ch);
                    }
                }

                // split
                {
                    int index = c.getIndex();
                    c.updateSplitFlag();
                    if (index > config.indelEndSkip && m - index > config.indelEndSkip) {
                        SearchState nextState = c.nextStateAfterSplit();
                        if (nextState != null)
                            queue.add(nextState);
                    }
                }

                assert (c.isFinished());
            }
        }

        private void reportAlignment(SearchState c) throws Exception {
            out.emit(c);
        }

        private void step(SearchState c, ACGT ch) {
            c.updateFlag(ch);
            Cursor nextCursor = c.cursor.next(ch);
            if (nextCursor == null)
                return; // no match

            SearchState nextState = c.nextState(ch, nextCursor, queryMask[c.getStrandIndex()]);
            if (nextState != null)
                queue.add(nextState);
        }

        private SearchState exactMatch(SearchState c) {
            Cursor cursor = c.cursor;
            final int n = cursor.getRemainingBases();
            final int strandIndex = c.getStrandIndex();
            int numExtend = 0;
            while (numExtend < n) {
                cursor = cursor.next(cursor.nextACGT(q));
                if (cursor == null)
                    return null;
                numExtend++;
            }
            return new SearchState(cursor, null, c.getNumDifferences());
        }

    }

}
