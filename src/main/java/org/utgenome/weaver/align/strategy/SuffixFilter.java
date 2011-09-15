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
// Since: 2011/09/15
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.strategy;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.PriorityQueue;

import org.utgenome.weaver.align.ACGT;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.AlignmentConfig;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.QueryMask;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.SuffixInterval;
import org.utgenome.weaver.align.record.AlignmentRecord;
import org.utgenome.weaver.align.record.Read;
import org.utgenome.weaver.align.record.SingleEndRead;
import org.utgenome.weaver.align.strategy.BidirectionalSuffixFilter.SFKey;
import org.utgenome.weaver.align.strategy.ReadAlignmentNFA.NextState;
import org.utgenome.weaver.parallel.Reporter;
import org.xerial.lens.SilkLens;
import org.xerial.util.log.Logger;

public class SuffixFilter
{
    private static Logger                   _logger               = Logger.getLogger(SuffixFilter.class);

    private final FMIndexOnGenome           fmIndex;
    private final AlignmentConfig           config;
    private final ACGTSequence              reference;

    /**
     * query length -> staircase filter of this query length
     */
    private HashMap<SFKey, StaircaseFilter> staircaseFilterHolder = new HashMap<SFKey, StaircaseFilter>();

    /**
     * Prepare a suffix filter
     * 
     * @param fmIndex
     *            FM index
     * @param config
     *            alignment score config
     * @param m
     *            read length
     */
    public SuffixFilter(FMIndexOnGenome fmIndex, ACGTSequence reference, AlignmentConfig config) {
        this.fmIndex = fmIndex;
        this.reference = reference;
        this.config = config;
    }

    public List<AlignmentRecord> align(ACGTSequence seq) throws Exception {
        return align(new SingleEndRead("read", seq, null));
    }

    public List<AlignmentRecord> align(Read read) throws Exception {
        final List<AlignmentRecord> alignmentResult = new ArrayList<AlignmentRecord>();
        new AlignmentProcess(read, new Reporter() {
            @Override
            public void emit(Object result) throws Exception {
                if (AlignmentRecord.class.isInstance(result)) {
                    AlignmentRecord r = (AlignmentRecord) result;
                    _logger.debug(SilkLens.toSilk("alignment", r));
                    alignmentResult.add(r);
                }
            }
        }).align();
        return alignmentResult;
    }

    public void align(Read read, Reporter out) throws Exception {
        new AlignmentProcess(read, out).align();
    }

    private class AlignmentProcess
    {
        private final Read     read;                                // original read
        private final int      m;                                   // read length
        private ACGTSequence[] q              = new ACGTSequence[2]; // forward/reverse query sequence
        private QueryMask[]    queryMask      = new QueryMask[2];   // bit mask of ACGT occurrence positions
        private Reporter       out;

        private final int      k;
        private int            minMismatches;
        private int            maxMatchLength = 0;
        private int            bestScore      = -1;

        private StateQueue     queue          = new StateQueue();

        public AlignmentProcess(Read read, Reporter out) {
            this.read = read;
            this.m = (int) read.getRead(0).textSize();
            this.q[0] = read.getRead(0);
            this.q[1] = q[0].complement();
            this.out = out;

            this.k = config.getMaximumEditDistances(m);
            this.minMismatches = k + 1;
        }

        /**
         * Get the staircase filter of NFA for the specified query length
         * 
         * @param queryLength
         * @return
         */
        private StaircaseFilter getStairCaseFilter(int queryLength) {
            int kk = minMismatches;
            SFKey key = new SFKey(kk, queryLength);
            if (!staircaseFilterHolder.containsKey(key)) {
                StaircaseFilter filter = new StaircaseFilter(queryLength, kk);
                staircaseFilterHolder.put(key, filter);
            }

            return staircaseFilterHolder.get(key);
        }

        private void initQueue(PrefixScan ps) {
            StaircaseFilter filter = getStairCaseFilter(m);
            final int s = filter.getNumChunks();
            for (int i = 0; i < s; ++i) {
                if (ps.chunkWithMismatch.get(i))
                    continue; // has mismatch

                int cs = filter.getChunkStart(i);
                int w = filter.getChunkSize(i);
                ReadAlignmentNFA nfa = new ReadAlignmentNFA(k);
                nfa.activateDiagonalStates();
                queue.add(new SFState(ps.strand, cs, cs + w, config.matchScore * w, ps.si.get(i), nfa, false));
            }

        }

        private int numFMIndexSearches = 0;
        private int numCutOff          = 0;
        private int numFiltered        = 0;

        public void align() {

            // Check whether the read contains too many Ns
            {
                long countN = q[0].fastCount(ACGT.N, 0, m);
                if (countN > k) {
                    return; // skip this alignment
                }

                if (countN > 0) {
                    q[0] = q[0].replaceN_withA();
                    q[1] = q[1].replaceN_withA();
                }
            }

            PrefixScan psF = PrefixScan.scanRead(fmIndex, q[0], Strand.FORWARD, getStairCaseFilter(m));
            PrefixScan psR = PrefixScan.scanRead(fmIndex, q[1], Strand.REVERSE, getStairCaseFilter(m));
            if (_logger.isDebugEnabled())
                _logger.debug("prefix scan: %s\t%s", psF, psR);

            initQueue(psF);
            initQueue(psR);

            queryMask[0] = new QueryMask(q[0]);
            queryMask[1] = new QueryMask(q[1]);

            StaircaseFilter filter = getStairCaseFilter(m);

            while (!queue.isEmpty()) {
                SFState c = queue.poll();

                if (_logger.isDebugEnabled()) {
                    _logger.debug("state: %s, #states:%d, FMSearch:%d, CutOff:%d, Filtered:%d", c, queue.size(),
                            numFMIndexSearches, numCutOff, numFiltered);
                }

                if (c.hasHit || c.index >= m || c.si.isUniqueHit()) {
                    addCandidate(c);
                    continue;
                }
                int nm = c.nfa.kOffset;
                if (c.scoreUpperBound(m) < bestScore || nm > minMismatches) {
                    numCutOff++;
                    continue;
                }

                final int strandIndex = c.strand.index;
                SuffixInterval[] nextSi = fmIndex.forwardSearch(c.strand, c.si);
                numFMIndexSearches++;
                for (ACGT ch : ACGT.exceptN) {
                    SuffixInterval si = nextSi[ch.code];
                    if (si != null) {
                        SFState nextState = c.nextState(ch, m, queryMask[strandIndex], si, getStairCaseFilter(m));
                        if (nextState != null) {
                            if (nextState.hasHit)
                                addCandidate(nextState);
                            else
                                queue.add(nextState);
                        }
                    }
                }
            }

        }

        public void addCandidate(SFState c) {
            _logger.debug("candidate:%s", c);
            if (c.si.isUniqueHit()) {

            }
            else {

            }
        }

    }

    private static class StateQueue extends PriorityQueue<SFState>
    {
        /**
         * 
         */
        private static final long serialVersionUID = 1L;

        public StateQueue() {
            super();
        }
    }

    private class SFState implements Comparable<SFState>
    {
        // |--------(read)-------|
        // |--(offset)|
        //               | <- index [0, m)
        public Strand                 strand;
        public final int              offset;
        public final int              index;
        public final int              score;
        public final SuffixInterval   si;
        public final ReadAlignmentNFA nfa;
        public final boolean          hasHit;

        public SFState(Strand strand, int offset, int index, int score, SuffixInterval si, ReadAlignmentNFA nfa,
                boolean hasHit) {
            this.strand = strand;
            this.offset = offset;
            this.index = index;
            this.score = score;
            this.si = si;
            this.nfa = nfa;
            this.hasHit = hasHit;
        }

        @Override
        public String toString() {
            return String.format("k%d%s%d/%d score:%s, si:%s", nfa.kOffset, strand.symbol, index, offset, score, si);
        }

        public int scoreUpperBound(int queryLen) {
            return score + (offset + (queryLen - index)) * config.matchScore;
        }

        public SFState nextState(ACGT ch, int queryLen, QueryMask queryMask, SuffixInterval nextSi,
                StaircaseFilter filter) {
            int nextIndex = index + 1;
            NextState next = nfa.nextState(nextIndex, queryLen, ch, queryMask, filter);
            if (next == null)
                return null;

            int diff = next.nextState.kOffset - nfa.kOffset;
            int newScore = score + score * config.mismatchPenalty;
            return new SFState(strand, offset, index + 1, newScore, nextSi, next.nextState, next.hasMatch);
        }

        @Override
        public int compareTo(SFState other) {
            return this.score - other.score;
        }
    }

}
