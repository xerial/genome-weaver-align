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
import org.utgenome.weaver.align.BitVector;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.Strand;
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
    private static Logger         _logger = Logger.getLogger(SuffixFilter.class);

    private final FMIndexOnGenome fmIndex;
    private final int             k;
    private final int             m;

    private StaircaseFilter       staircaseFilter;

    /**
     * NFA
     * 
     * @author leo
     * 
     */
    public static class SearchState
    {
        private final BitVector[] automaton;
        // 32 bit = strand(1) + searchDirection(2) + searchFlag (5) + minK (8) + index (16)  
        private int               state = 0;

        private SearchState(BitVector[] automaton, Strand strand, SearchDirection searchDirection, int minK, int index) {
            this.automaton = automaton;
            this.state = (strand.index & 1) | ((searchDirection.index & 0x03) << 1);
            this.state = newFlag(minK, index);
        }

        private SearchState(BitVector[] automaton, int newState) {
            this.automaton = automaton;
            this.state = newState;
        }

        public int getIndex() {
            return (state >>> 16) & 0xFFFF;
        }

        public int getMinDifferences() {
            return (state >>> 8) & 0xFF;
        }

        public boolean isFinished() {
            return ((state >>> 3) & 0x1F) == 0x1F; // ACGT + split
        }

        private int newFlag(int minDiff, int newIndex) {
            return (state & 0x07) | ((minDiff & 0xFF) << 8) | ((newIndex & 0xFFFF) << 16);
        }

        public static SearchState initialState(Strand strand, SearchDirection searchDirection, int k, int m) {
            BitVector[] automaton = new BitVector[k + 1];
            // Activate the diagonal states 
            for (int i = 0; i <= k; ++i) {
                automaton[i] = new BitVector(m);
                // TODO optimize flag set
                for (int j = 0; j <= i; ++j)
                    automaton[i].set(j);
            }
            SearchState s = new SearchState(automaton, strand, searchDirection, 0, 0);
            return s;
        }

        public SearchState nextStateAfterSplit(StaircaseFilter staircaseFilter) {
            this.state |= 1 << 5;
            // use the same automaton state
            final int k = automaton.length - 1;
            if (getMinDifferences() < k) {
                return new SearchState(automaton, newFlag(getMinDifferences() + 1, getIndex()));
            }
            else
                return null;
        }

        public SearchState nextState(ACGT ch, QueryMask queryMask, StaircaseFilter staircaseFilter) {
            // update the search flag
            this.state |= 1 << ch.code;

            final int k = automaton.length - 1;
            final int m = (int) automaton[0].size();

            BitVector[] prev = automaton;
            BitVector[] next = new BitVector[k + 1];

            final BitVector qeq = queryMask.getPatternMask(ch);

            final int nextIndex = getIndex() + 1;

            int nm = 0;
            // R'_0 = (R_0 << 1) | P[ch]
            next[nm] = prev[nm].rshift(1)._and(qeq);
            if (next[nm].get(m - 1)) {
                // Found a full match
                //return new DFSState(State.Finished, 0);
                return new SearchState(next, newFlag(nm, nextIndex));
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
                    return new SearchState(next, newFlag(i, nextIndex));

                if (nm == i && !next[i].get(nextIndex))
                    ++nm;
            }

            if (nm >= k) {
                // no match
                return null;
            }
            else {
                // extend the match
                return new SearchState(next, newFlag(nm, nextIndex));
            }
        }
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

    public SuffixFilter(FMIndexOnGenome fmIndex, int k, int m) {
        this.fmIndex = fmIndex;
        this.k = k;
        this.m = m;
        this.staircaseFilter = new StaircaseFilter(m, k);
    }

    public void align(ACGTSequence query, Reporter out) {
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

        public void align(Reporter out) {
            // Add states for both strands
            queue.add(SearchState.initialState(Strand.FORWARD, SearchDirection.Forward, k, m));
            queue.add(SearchState.initialState(Strand.REVERSE, SearchDirection.Forward, k, m));

            while (!queue.isEmpty()) {

            }
        }

    }

}
