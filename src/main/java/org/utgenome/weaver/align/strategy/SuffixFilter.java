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
     * Cursor for bidirectional search
     * 
     * @author leo
     * 
     */
    public class Cursor
    {
        // flag(8bit) :=  strand(1), searchDirection(2), extension type(2)
        private final byte flag;
        private Score      score;
        public final int   cursorF;
        public final int   cursorB;

        Cursor(byte flag, Score score, int cursorF, int cursorB) {
            this.flag = flag;
            this.score = score;
            this.cursorF = cursorF;
            this.cursorB = cursorB;
        }

        public Cursor(Strand strand, SearchDirection searchDirection, Score score, ExtensionType extensionType,
                int cursorF, int cursorB) {
            this((byte) (strand.index | (searchDirection.index << 1) | (extensionType.code << 3)), score, cursorF,
                    cursorB);
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

            // increment num split
            return new Cursor(flag, score.extendWithSplit(config), cursor, cursor);
        }

        public Cursor next() {
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
                if (nextF < m) {
                    ++nextF;
                }
                else {
                    // switch to backward search
                    d = SearchDirection.Backward;
                    --nextB;
                }
                break;
            }

            return new Cursor(flag, score, nextF, nextB);
        }
    }

    /**
     * NFA
     * 
     * @author leo
     * 
     */
    public class SearchState
    {
        public final Cursor         cursor;
        private final BitVector[]   automaton;
        // 32 bit = searchFlag (5) + minK (8)
        private int                 state;

        public final SuffixInterval siF;
        public final SuffixInterval siB;

        private SearchState(Cursor cursor, BitVector[] automaton, int minK, SuffixInterval siF, SuffixInterval siB) {
            this.cursor = cursor;
            this.automaton = automaton;
            this.state = minK << 5;
            this.siF = siF;
            this.siB = siB;
        }

        public int getStrandIndex() {
            return cursor.flag & 1;
        }

        public int getNumDifferences() {
            return this.cursor.score.layer();
        }

        public int getIndex() {
            return cursor.getIndex();
        }

        public int getMinDifferences() {
            return (state >>> 5) & 0xFF;
        }

        public boolean isFinished() {
            return (state & 0x1F) == 0x1F; // ACGT + split
        }

        public void updateFlag(ACGT ch) {
            this.state |= 1 << ch.code;
        }

        public void updateSplitFlag() {
            this.state |= 1 << 5;
        }

        public boolean isChecked(ACGT ch) {
            return (state & (1 << ch.code)) != 0;
        }

        public SearchState nextStateAfterSplit() {
            updateSplitFlag();
            // use the same automaton state
            if (getMinDifferences() < k) {
                return new SearchState(cursor.split(), automaton, getMinDifferences() + 1, fmIndex.wholeSARange(),
                        fmIndex.wholeSARange());
            }
            else
                return null;
        }

        public SearchState nextState(ACGT ch, SuffixInterval siF, SuffixInterval siB, QueryMask queryMask) {
            // update the search flag
            updateFlag(ch);

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
                return new SearchState(cursor.next(), next, nm, siF, siB);
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
                    return new SearchState(cursor.next(), next, i, siF, siB);

                if (nm == i && !next[i].get(nextIndex))
                    ++nm;
            }

            if (nm >= k) {
                // no match
                return null;
            }
            else {
                // extend the match
                return new SearchState(cursor.next(), next, nm, siF, siB);
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
        SearchState s = new SearchState(
                new Cursor(strand, searchDirection, Score.initial(), ExtensionType.MATCH, 0, 0), automaton, 0,
                fmIndex.wholeSARange(), fmIndex.wholeSARange());
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
                patternMask[i] = new BitVector(m);

            for (int i = 0; i < m; ++i) {
                int index = offset + i;
                if (index >= m) {
                    // for bidirectional search
                    index = m - i - 1;
                }
                ACGT ch = query.getACGT(index);
                if (ch == ACGT.N) {
                    for (ACGT each : ACGT.exceptN)
                        patternMask[ch.code].set(i);
                }
                else
                    patternMask[ch.code].set(i);
            }
        }

        public BitVector getPatternMask(ACGT ch) {
            return patternMask[ch.code];
        }
    }

    public SuffixFilter(FMIndexOnGenome fmIndex, AlignmentScoreConfig config, int m) {
        this.fmIndex = fmIndex;
        this.config = config;
        this.k = config.maximumEditDistances;
        this.m = m;
        this.staircaseFilter = new StaircaseFilter(m, k);
    }

    public void align(ACGTSequence query, Reporter out) throws Exception {
        new AlignmentProcess(query).align(out);
    }

    private static class StatePreference implements Comparator<SearchState>
    {
        @Override
        public int compare(SearchState o1, SearchState o2) {
            // prefer longer match
            int diff = o2.getIndex() - o1.getIndex();
            if (diff != 0)
                return diff;

            int kDiff = o1.getMinDifferences() - o2.getMinDifferences();
            return kDiff;
        }
    }

    public class AlignmentProcess
    {

        private ACGTSequence[]     q         = new ACGTSequence[2];
        private QueryMask[]        queryMask = new QueryMask[2];
        PriorityQueue<SearchState> queue     = new PriorityQueue<SearchState>(11, new StatePreference());

        public AlignmentProcess(ACGTSequence query) {
            q[0] = query;
            q[1] = query.complement();
            queryMask[0] = new QueryMask(q[0]);
            queryMask[1] = new QueryMask(q[1]);
        }

        public void align(Reporter out) throws Exception {
            // Add states for both strands
            queue.add(initialState(Strand.FORWARD, SearchDirection.Forward, k, m));
            queue.add(initialState(Strand.REVERSE, SearchDirection.Forward, k, m));

            while (!queue.isEmpty()) {
                SearchState c = queue.peek();
                if (c.isFinished()) {
                    queue.poll();
                    continue;
                }

                int strandIndex = c.getStrandIndex();
                int allowedMismatches = k - c.getNumDifferences();
                if (c.getIndex() + allowedMismatches >= m) {
                    // found a hit
                    out.emit(c);
                    queue.poll();
                    continue;
                }

                if (allowedMismatches == 0) {
                    // exact match
                }

                ACGT nextBase = c.cursor.nextACGT(q);
                if (!c.isChecked(nextBase)) {
                    // search for a base in the read
                    c.updateFlag(nextBase);
                    BidirectionalSuffixInterval nextSi = extendSearch(c, nextBase);
                    if (nextSi != null) {
                        enqueue(c.nextState(nextBase, nextSi.forwardSi, nextSi.backwardSi, queryMask[strandIndex]));
                    }
                    continue;
                }

                // search for the other bases
                for (ACGT ch : ACGT.exceptN) {
                    if (ch == nextBase)
                        continue;

                    if (!c.isChecked(ch)) {
                        c.updateFlag(ch);
                        BidirectionalSuffixInterval nextSi = extendSearch(c, ch);
                        enqueue(c.nextState(ch, nextSi.forwardSi, nextSi.backwardSi, queryMask[strandIndex]));
                    }
                }

                // split
                {
                    c.updateSplitFlag();
                    enqueue(c.nextStateAfterSplit());
                }

            }
        }

        private void enqueue(SearchState c) {
            if (c == null)
                return;

            queue.add(c);
        }

        BidirectionalSuffixInterval extendSearch(SearchState c, ACGT ch) {

            SuffixInterval siF = c.siF;
            SuffixInterval siB = c.siB;

            Strand strand = c.cursor.getStrand();

            switch (c.cursor.getSearchDirection()) {
            case Forward:
                siF = fmIndex.forwardSearch(strand, ch, siF);
                if (siF.isEmpty())
                    return null;
                break;
            case Backward:
                siB = fmIndex.backwardSearch(strand, ch, siB);
                if (siB.isEmpty())
                    return null;
                break;
            case BidirectionalForward:
                if (c.cursor.cursorF < m) {
                    BidirectionalSuffixInterval bSi = fmIndex.bidirectionalForwardSearch(strand, ch,
                            new BidirectionalSuffixInterval(siF, siB));
                    if (bSi == null)
                        return null;
                    siF = bSi.forwardSi;
                    siB = bSi.backwardSi;
                }
                else {
                    siB = fmIndex.backwardSearch(strand, ch, siB);
                    if (siB == null)
                        return null;
                }
                break;
            }

            return new BidirectionalSuffixInterval(siF, siB);

        }

    }

}
