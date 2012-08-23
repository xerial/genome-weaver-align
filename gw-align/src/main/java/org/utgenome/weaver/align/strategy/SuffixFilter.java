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
import java.util.List;
import java.util.PriorityQueue;
import java.util.TreeSet;

import org.utgenome.UTGBException;
import org.utgenome.weaver.align.ACGT;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.Aligner;
import org.utgenome.weaver.align.AlignmentConfig;
import org.utgenome.weaver.align.BitParallelSmithWaterman;
import org.utgenome.weaver.align.CIGAR;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.QueryMask;
import org.utgenome.weaver.align.SequenceBoundary.PosOnGenome;
import org.utgenome.weaver.align.SmithWatermanAligner.Alignment;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.SuffixInterval;
import org.utgenome.weaver.align.record.AlignmentRecord;
import org.utgenome.weaver.align.record.Read;
import org.utgenome.weaver.align.record.ReadHit;
import org.utgenome.weaver.align.record.SingleEndRead;
import org.utgenome.weaver.align.strategy.ReadAlignmentNFA.NextState;
import org.utgenome.weaver.align.strategy.StaircaseFilter.StaircaseFilterHolder;
import org.utgenome.weaver.parallel.Reporter;
import org.xerial.lens.SilkLens;
import org.xerial.util.StopWatch;
import org.xerial.util.log.Logger;

public class SuffixFilter implements Aligner
{
    private static Logger         _logger               = Logger.getLogger(SuffixFilter.class);

    private final FMIndexOnGenome fmIndex;
    private final AlignmentConfig config;
    private final ACGTSequence    reference;

    /**
     * query length -> staircase filter of this query length
     */
    private StaircaseFilterHolder staircaseFilterHolder = StaircaseFilter.newHolder();

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

            this.k = config.getMaximumEditDistance(m);
            this.minMismatches = k + 1;
        }

        /**
         * Get the staircase filter of NFA for the specified query length
         * 
         * @param queryLength
         * @return
         */
        private StaircaseFilter getStairCaseFilter(int queryLength) {
            return staircaseFilterHolder.getStairCaseFilter(queryLength, minMismatches);
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
        private int numSW              = 0;

        public void align() throws Exception {

            StopWatch s = new StopWatch();
            try {
                align_internal();
                boolean hasHit = minMismatches <= k && !resultHolder.hitList.isEmpty();
                ReadHit besthit = null;
                String cigar = "";
                if (hasHit) {
                    besthit = resultHolder.hitList.get(0);
                    cigar = besthit.getCigarConcatenated();
                }
                boolean isUnique = resultHolder.isUnique();

                if (_logger.isDebugEnabled()) {
                    _logger.debug("query:%s %s %2s %s k:%d, FM Search:%,d, SW:%d, CutOff:%d, Filtered:%d, %.5f sec.",
                            read.name(), hasHit ? besthit.strand.symbol : " ", hasHit ? besthit.getAlignmentState()
                                    : " ", cigar, minMismatches, numFMIndexSearches, numSW, numCutOff, numFiltered, s
                                    .getElapsedTime());

                    if (minMismatches > k || numFMIndexSearches > 500) {
                        _logger.debug("query:%s", q[0]);
                    }
                }
                if (_logger.isTraceEnabled())
                    _logger.trace("qual :%s", read.getQual(0));

                // Issue 28 
                if (hasHit) {
                    switch (config.reportType) {
                    case ALLHITS:
                        for (ReadHit each : resultHolder.hitList) {
                            report(each, 0);
                        }
                        break;
                    case BESTHIT:
                        report(besthit, 0);
                        break;
                    case TOPL: {
                        int max = Math.min(config.topL, resultHolder.hitList.size());
                        for (int i = 0; i < max; ++i) {
                            report(resultHolder.hitList.get(i), 0);
                        }
                        break;
                    }
                    }
                }
                else {
                    // report unmapped read
                    report(new ReadHit("*", 0, 0, 0, 0, -1, Strand.FORWARD, new CIGAR(), 0, null), 0);
                }

            }
            catch (Exception e) {
                _logger.error("error at query: %s", q[0]);
                throw e;
            }
        }

        public void report(ReadHit hit, int numTotalHits) throws Exception {
            AlignmentRecord r = AlignmentRecord.convert(hit, read, numTotalHits);
            out.emit(r);
        }

        public void align_internal() {

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
            if (_logger.isTraceEnabled())
                _logger.debug("prefix scan: %s\t%s", psF, psR);

            initQueue(psF);
            initQueue(psR);

            queryMask[0] = new QueryMask(q[0]);
            queryMask[1] = new QueryMask(q[1]);

            StaircaseFilter filter = getStairCaseFilter(m);

            while (!queue.isEmpty()) {
                SFState c = queue.poll();

                int nm = c.nfa.kOffset;
                if (nm > minMismatches || c.scoreUpperBound(m) < bestScore) {
                    numCutOff++;
                    continue;
                }

                if (_logger.isTraceEnabled()) {
                    _logger.trace("state: %s, #states:%d, FMSearch:%d, SW:%d, CutOff:%d, Filtered:%d", c, queue.size(),
                            numFMIndexSearches, numSW, numCutOff, numFiltered);
                }

                if (c.hasHit || c.index >= m || c.si.isUniqueHit()) {
                    addCandidate(c);
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
                            queue.add(nextState);
                        }
                        else
                            numFiltered++;
                    }
                }
            }

        }

        @SuppressWarnings("unchecked")
        private TreeSet<Long>[] candidates = new TreeSet[] { new TreeSet<Long>(), new TreeSet<Long>() };
        private List<ReadHit>   hitList    = new ArrayList<ReadHit>();

        void addCandidate(SFState c) {

            if (c.si.isUniqueHit()) {
                long seqIndex = fmIndex.toCoordinate(c.si.lowerBound, c.strand, SearchDirection.Forward);
                int offsetFromSearchHead = c.strand.isForward() ? c.index : m - c.index;
                long start = seqIndex - offsetFromSearchHead;

                //                if (_logger.isTraceEnabled())
                //                    _logger.trace("candidate seq index:%d, pos:%d", seqIndex, start);

                if (candidates[c.strand.index].contains(start))
                    return;
                else
                    candidates[c.strand.index].add(start);

                // verify the newly added entry
                long refStart = Math.max(0, start - k);
                long refEnd = Math.min(start + m + k, fmIndex.textSize());

                ACGTSequence ref = reference.subString(refStart, refEnd);
                ACGTSequence query = q[c.strand.index];
                if (!c.strand.isForward())
                    query = query.reverse();
                Alignment alignment = BitParallelSmithWaterman.alignBlockDetailed(ref, query, config.bandWidth);
                numSW++;
                if (alignment != null) {

                    if (_logger.isDebugEnabled() && alignment.numMismatches < minMismatches)
                        _logger.debug("Found an alignment %s", alignment);

                    try {
                        PosOnGenome p = fmIndex.getSequenceBoundary().translate(refStart + alignment.pos + 1,
                                Strand.FORWARD);
                        reportResult(new ReadHit(p.chr, p.pos, m, 0, m, alignment.numMismatches, c.strand,
                                alignment.cigar, (int) c.si.range(), null));
                    }
                    catch (UTGBException e) {
                        _logger.error(e);
                    }
                }
            }
            else {
                // multi hit

            }
        }

        void reportResult(ReadHit hit) {
            if (hit.getTotalMatchLength() == 0)
                return; // no match

            if (hit.diff > minMismatches)
                return;

            resultHolder.add(hit);
        }

        private AlignmentResultHolder resultHolder = new AlignmentResultHolder();

        private class AlignmentResultHolder
        {
            ArrayList<ReadHit> hitList = new ArrayList<ReadHit>();

            public int totalHits() {
                int count = 0;
                for (ReadHit each : hitList) {
                    count += each.numHits;
                }
                return count;
            }

            public void add(ReadHit hit) {
                int newK = hit.getTotalDifferences();
                int matchLen = hit.getTotalMatchLength();
                if (matchLen > 0 && newK < minMismatches) {
                    minMismatches = newK;
                    _logger.trace("min mismatches: %d", minMismatches);
                    bestScore = hit.getTotalScore(config);
                }
                if (maxMatchLength < matchLen) {
                    maxMatchLength = matchLen;
                }

                ArrayList<ReadHit> newList = new ArrayList<ReadHit>();
                for (ReadHit each : hitList) {
                    if (each.getTotalDifferences() <= minMismatches) {
                        newList.add(each);
                    }
                }
                newList.add(hit);
                hitList = newList;
            }

            public boolean isUnique() {
                if (hitList.size() == 1) {
                    return hitList.get(0).numHits == 1;
                }
                return false;
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
            return String.format("%sk%d%s%d/%d score:%d, si:%s", hasHit ? "*" : " ", nfa.kOffset, strand.symbol, index,
                    offset, score, si);
        }

        public int scoreUpperBound(int queryLen) {
            return score + (offset + (queryLen - index)) * config.matchScore;
        }

        public SFState nextState(ACGT ch, int queryLen, QueryMask queryMask, SuffixInterval nextSi,
                StaircaseFilter filter) {
            int nextIndex = index + 1;
            NextState next = nfa.nextState(nextIndex, nextIndex - offset, queryLen - offset, ch, queryMask, filter);
            if (next == null)
                return null;

            int diff = next.nextState.kOffset - nfa.kOffset;
            int newScore = score - diff * config.mismatchPenalty;
            if (diff == 0)
                newScore++;
            return new SFState(strand, offset, nextIndex, newScore, nextSi, next.nextState, next.hasMatch);
        }

        @Override
        public int compareTo(SFState other) {
            int diff = (nfa.kOffset - other.nfa.kOffset);
            if (diff == 0)
                diff = -(this.score - other.score);

            return diff;
        }
    }

}
