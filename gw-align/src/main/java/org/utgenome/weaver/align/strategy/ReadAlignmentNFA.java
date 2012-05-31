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
// FMSearchNFA.java
// Since: 2011/08/25
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.strategy;

import org.utgenome.weaver.align.ACGT;
import org.utgenome.weaver.align.QueryMask;

/**
 * Non-deterministic finite automaton holding search states of FM-index.
 * 
 * <pre>
 *     q0 q1 q2 ... q_{m-1}
 *  k=0 +--+--+--
 *      |\ |\ |\
 *      | \| \| \  ...
 *  k=1 +--+--+--
 * 
 * </pre>
 * 
 * The states in the automaton are represented as bit flags. kOffset is used to
 * reduce the size of automaton. When no states become active at some layer k,
 * the layer will not be used in future, so we can remove that layer from the
 * automaton.
 * 
 * <h3>transitions</h3>
 * <ul>
 * <li>right:match
 * <li>down:deletion from reference
 * <li>right down: insertion to reference (epsilon transition) or mismatch
 * </ul>
 * 
 * @author leo
 * 
 */
public class ReadAlignmentNFA
{
    // bit flags holding states at column [index - (k-kOffset), index + (k-kOffset) + 1], where k is # of allowed mismatches 
    private final long[] automaton;
    public final int     kOffset;

    public ReadAlignmentNFA(int numAllowedMismatches) {
        this(new long[numAllowedMismatches + 1], 0);
    }

    private ReadAlignmentNFA(long[] automaton, int kOffset) {
        this.automaton = automaton;
        this.kOffset = kOffset;
    }

    @Override
    public String toString() {
        return toNFAStateString();
    }

    public ReadAlignmentNFA activateDiagonalStates() {
        // Activate the diagonal states 
        final int k = automaton.length - 1;
        for (int i = 0; i < automaton.length; ++i) {
            automaton[i] = 1L << (k + i);
        }
        return this;
    }

    private String toBinary(long val, int w) {
        StringBuilder s = new StringBuilder();
        for (int i = 0; i < w; ++i) {
            s.append(((val >>> i) & 1L) == 0 ? "0" : "1");
        }
        return s.toString();
    }

    public String toNFAStateString() {
        int w = (this.automaton.length + 1) * 2 + 1;
        StringBuilder s = new StringBuilder();
        s.append(String.format("kOffset:%d\n", kOffset));
        for (int j = 0; j < automaton.length; ++j) {
            s.append(toBinary(automaton[j], w));
            if (j != automaton.length - 1) {
                s.append(" \n");
            }
        }
        return s.toString();
    }

    public ReadAlignmentNFA nextStateAfterSplit(int k) {
        int height = automaton.length - 1;
        long[] nextAutomaton = new long[height];
        for (int i = 0; i < height; ++i) {
            nextAutomaton[i] = 1L << (height + i - 1);
        }
        return new ReadAlignmentNFA(nextAutomaton, kOffset + 1);
    }

    public static class NextState
    {
        public final ReadAlignmentNFA nextState;
        public final boolean          hasMatch;

        public NextState(ReadAlignmentNFA nextState, boolean hasMatch) {
            this.nextState = nextState;
            this.hasMatch = hasMatch;
        }
    }

    public NextState nextState(Cursor cursor, ACGT ch, QueryMask queryMask, StaircaseFilter staircaseFilter) {
        final int height = automaton.length;
        final int k = kOffset + height - 1;
        final int kr = k - kOffset;
        final long qeq = queryMask.getBidirectionalPatternMask64(cursor.getSearchDirection(),
                cursor.getNextACGTIndex(), cursor.pivot, cursor.cursor, ch, kr);
        return nextState(qeq, cursor.getProcessedBases(), cursor.getFragmentLength(), ch, queryMask, staircaseFilter);
    }

    public NextState nextState(int nextACGTIndex, int progressIndex, int m, ACGT ch, QueryMask queryMask,
            StaircaseFilter staircaseFilter) {
        final int height = automaton.length;
        final int k = kOffset + height - 1;
        final int kr = k - kOffset;
        final long qeq = queryMask.getBidirectionalPatternMask64(SearchDirection.Forward, nextACGTIndex, 0,
                nextACGTIndex, ch, kr);
        return nextState(qeq, progressIndex, m, ch, queryMask, staircaseFilter);
    }

    private NextState nextState(long qeq, int progressIndex, int fragmentLength, ACGT ch, QueryMask queryMask,
            StaircaseFilter staircaseFilter) {
        final int height = automaton.length;
        final int k = kOffset + height - 1;
        final int kr = k - kOffset;

        long[] prev = automaton;
        long[] next = new long[automaton.length];

        int minKwithMatch = k + 1;
        int minKwithProgress = k + 1;

        String qeqStr = toBinary(qeq, 10);

        //List<Integer> matchPos = new ArrayList<Integer>();
        // Update the automaton
        // R'_0 = ((R_0 & P[ch]) << 1) & (suffix filter)
        next[0] = (prev[0] & qeq) << 1;
        if (next[0] != 0) {
            minKwithMatch = 0;
            minKwithProgress = 0;
            //matchPos.add(cursor.getNextACGTIndex());
        }

        next[0] &= staircaseFilter.getStairCaseMask64bit(kOffset, progressIndex - k);
        for (int i = 1; i < height; ++i) {
            // R'_{i+1} = ((R_{i+1} & P[ch]) << 1) | R_i | (R_i << 1) | (R'_i << 1) 
            next[i] = (prev[i] & qeq) << 1;
            if (minKwithMatch > k && next[i] != 0) {
                minKwithMatch = i;

            }
            next[i] |= prev[i - 1] | (prev[i - 1] << 1) | (next[i - 1] << 1);
            // Apply a suffix filter (staircase mask)
            next[i] &= staircaseFilter.getStairCaseMask64bit(kOffset + i, progressIndex - k);
            if (minKwithProgress > k && (next[i] & (1L << height)) != 0L) {
                minKwithProgress = i;
            }
        }

        // Find a match at query position m
        final int m = fragmentLength;
        final int mPos = k + m - progressIndex;
        if (mPos < 64) {
            for (int nm = 0; nm < height; ++nm) {
                if ((next[nm] & (1L << mPos)) != 0L) {
                    return new NextState(new ReadAlignmentNFA(removeLayersFromAutomaton(next, nm), kOffset + nm), true);
                }
            }
        }

        // 
        int minK = Math.min(minKwithMatch, minKwithProgress);
        if (minK < k)
            return new NextState(new ReadAlignmentNFA(removeLayersFromAutomaton(next, minK), kOffset + minK), false);

        return null;
    }

    private long[] removeLayersFromAutomaton(long[] next, int numLayersToRemove) {
        if (numLayersToRemove == 0)
            return automaton;
        int nextHeight = automaton.length - numLayersToRemove;
        // trim rows with fewer mismatches  
        long[] trimmed = new long[nextHeight];
        for (int h = 0; h < nextHeight; ++h) {
            trimmed[h] = next[h + numLayersToRemove] >>> 1;
        }
        return trimmed;
    }

}
