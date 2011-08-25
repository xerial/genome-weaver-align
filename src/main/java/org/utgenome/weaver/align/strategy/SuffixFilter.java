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

import org.utgenome.gwt.utgb.client.bio.CIGAR;
import org.utgenome.weaver.align.ACGT;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.AlignmentScoreConfig;
import org.utgenome.weaver.align.BitParallelSmithWaterman;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.QueryMask;
import org.utgenome.weaver.align.SequenceBoundary.PosOnGenome;
import org.utgenome.weaver.align.SiSet;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.SuffixInterval;
import org.utgenome.weaver.align.record.AlignmentRecord;
import org.utgenome.weaver.align.record.Read;
import org.utgenome.weaver.align.record.ReadHit;
import org.utgenome.weaver.align.record.SWResult;
import org.utgenome.weaver.align.record.SingleEndRead;
import org.utgenome.weaver.align.strategy.FMSearchNFA.NextState;
import org.utgenome.weaver.parallel.Reporter;
import org.xerial.lens.SilkLens;
import org.xerial.util.StopWatch;
import org.xerial.util.StringUtil;
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
    private final ACGTSequence         reference;
    private final int                  k;
    private final int                  m;

    private StaircaseFilter            staircaseFilter;

    public static enum Diff {
        Match, Mismatch, Insertion, Deletion, Split
    };

    /**
     * NFA state and alignment cursor holder.
     * 
     * This NFA holds 2k+1 columns and k rows
     * 
     * @author leo
     * 
     */
    public class SearchState
    {
        public final Cursor         cursor;
        public final SuffixInterval currentSi;
        public final SiSet          siTable;
        private FMSearchNFA         automaton;
        // 32 bit = searchFlag(5) + currentBase(3) + minK (8) + priority(8) + hasHit(1) 
        private int                 state;
        public final SearchState    split;

        public int getLowerBoundOfK() {
            return (state >>> 8) & 0xFF;
        }

        public void setLowerBoundOfK(int diff) {
            state &= ~(0xFF << 8);
            state |= (diff & 0xFF) << 8;
        }

        /**
         * Lower value has higher priority
         * 
         * @return
         */
        public int getPriority() {
            return (state >>> 16) & 0xFF;
        }

        public void lowerThePrioity(int margin) {
            int priority = getPriority() + margin;
            state &= ~(0xFF << 16);
            state |= (priority & 0xFF) << 16;
        }

        public boolean hasHit() {
            return ((state >>> 24) & 1L) != 0;
        }

        public boolean isFinished() {
            return (state & 0x1F) == 0x1F; // ACGT + split
        }

        public ACGT currentACGT() {
            return ACGT.decode((byte) ((state >>> 5) & 0x7));
        }

        public void updateFlag(ACGT ch) {
            this.state |= 1 << ch.code;
        }

        public void fillSearchFlags() {
            this.state |= 0x1F;
        }

        public void updateSplitFlag() {
            this.state |= 1 << 4;
        }

        public boolean isSplitChecked() {
            return (state & (1 << 4)) != 0;
        }

        public boolean isChecked(ACGT ch) {
            return (state & (1 << ch.code)) != 0;
        }

        SearchState(SuffixInterval currentSi, ACGT ch, Cursor cursor, SiSet si, FMSearchNFA automaton,
                boolean hasMatch, int minK, int expectedK, SearchState split) {
            this.currentSi = currentSi;
            this.cursor = cursor;
            this.siTable = si;
            this.automaton = automaton;
            this.state = ((ch.code & 0x7) << 5) | ((minK & 0xFF) << 8) | ((expectedK & 0xFF) << 16)
                    | ((hasMatch ? 1 : 0) << 24);
            this.split = split;
        }

        /**
         * Initial forward search state
         * 
         * @param strand
         * @param searchDirection
         */
        public SearchState(SuffixInterval currentSi, Cursor cursor, int expectedK) {
            this(currentSi, ACGT.N, cursor, fmIndex.initSet(cursor.getSearchDirection()), new FMSearchNFA(k), false, 0,
                    expectedK, null);
            automaton.activateDiagonalStates();
        }

        String getUpdateFlag() {
            StringBuilder s = new StringBuilder();
            for (ACGT ch : ACGT.exceptN) {
                s.append(isChecked(ch) ? ch.name() : ch.name().toLowerCase());
            }
            if (isSplitChecked())
                s.append("^");
            return s.toString();
        }

        @Override
        public String toString() {
            return String.format("%sk%dp%d%s %s%s", hasHit() ? "*" : "", getLowerBoundOfK(), getPriority(), cursor,
                    getUpdateFlag(), siTable);
        }

        public int score() {
            int nm = getLowerBoundOfK();
            int mm = cursor.getIndex() - nm;
            return mm * config.matchScore - nm * config.mismatchPenalty;
        }

        public SearchState nextStateAfterSplit(ACGT ch) {
            updateSplitFlag();
            // use the same automaton state
            int minK = getLowerBoundOfK();
            if (minK < k) {
                return new SearchState(fmIndex.wholeSARange(), ch, cursor.split(), fmIndex.initSet(cursor
                        .getSearchDirection()), automaton.nextStateAfterSplit(k), false, minK, getPriority(), this);
            }
            else
                return null;
        }

        public SearchState nextState(ACGT ch, SiSet nextSi, QueryMask queryMask) {

            NextState next = automaton.nextState(cursor, ch, queryMask, staircaseFilter);
            if (next == null)
                return null; // no match

            int index = cursor.getIndex();
            int nextMinK = next.nextState.kOffset;

            SuffixInterval si = cursor.isForwardSearch() ? this.siTable.getForward(ch) : this.siTable.getBackward(ch);
            return new SearchState(si, ch, cursor.next(), nextSi, next.nextState, next.hasMatch, nextMinK,
                    getPriority(), split);
        }

    }

    private SiSet next(SearchState c, ACGT ch) {
        return c.cursor.nextSi(fmIndex, c.siTable, ch);
    }

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
    public SuffixFilter(FMIndexOnGenome fmIndex, ACGTSequence reference, AlignmentScoreConfig config, long m) {
        this.fmIndex = fmIndex;
        this.reference = reference;
        this.config = config;
        this.k = config.maximumEditDistances;
        this.m = (int) m;
        this.staircaseFilter = new StaircaseFilter(this.m, k);
    }

    public void align(ACGTSequence seq) throws Exception {
        align(new SingleEndRead("read", seq, null));
    }

    public void align(Read read) throws Exception {
        new AlignmentProcess(read, new Reporter() {
            @Override
            public void emit(Object result) throws Exception {
                _logger.debug(SilkLens.toSilk("result", result));

            }
        }).align();
    }

    public void align(Read read, Reporter out) throws Exception {
        new AlignmentProcess(read, out).align();
    }

    private static class StateQueue extends PriorityQueue<SearchState>
    {
        private static final long serialVersionUID = 1L;

        public StateQueue() {
            super(11, new Comparator<SearchState>() {
                // Comparator for selecting next target state to search
                @Override
                public int compare(SearchState o1, SearchState o2) {
                    int diff = 0;
                    diff = o1.getPriority() - o2.getPriority();
                    if (diff == 0)
                        diff = -(o1.score() - o2.score());
                    if (diff == 0)
                        diff = -(o1.cursor.getIndex() - o2.cursor.getIndex());

                    return diff;
                }
            });
        }

        @Override
        public String toString() {
            return StringUtil.join(this, "\n");
        }
    }

    /**
     * Alignment procedure
     * 
     * @author leo
     * 
     */
    class AlignmentProcess
    {

        private final Read     read;
        private ACGTSequence[] q             = new ACGTSequence[2];
        private QueryMask[]    queryMask     = new QueryMask[2];
        private StateQueue     queue         = new StateQueue();

        private Reporter       out;

        private int            minMismatches = k + 1;

        public AlignmentProcess(Read read, Reporter out) {
            this.read = read;
            this.q[0] = read.getRead(0);
            this.q[1] = q[0].complement();
            this.out = out;
        }

        public void align() throws Exception {

            StopWatch s = new StopWatch();
            try {
                align_internal();

                if (_logger.isDebugEnabled()) {
                    _logger.debug("query:%s - %s min K:%d, FM Search:%,d, Exact:%d, CutOff:%d, Filtered:%d, %.5f sec.",
                            read.name(), minMismatches <= k ? "(*)" : "   ", minMismatches, numFMIndexSearches,
                            numExactSearchCount, numCutOff, numFiltered, s.getElapsedTime());

                    if (numFMIndexSearches > 500) {
                        _logger.debug("query:%s", q[0]);
                    }
                }
                if (_logger.isTraceEnabled())
                    _logger.trace("qual :%s", read.getQual(0));

            }
            catch (Exception e) {
                _logger.error("error at query: %s", q[0]);
                throw e;
            }
        }

        /**
         * @throws Exception
         */
        /**
         * @throws Exception
         */
        public void align_internal() throws Exception {

            // Check whether the read contains too many Ns  
            if (q[0].fastCount(ACGT.N, 0, m) > config.maximumEditDistances) {
                return; // skip this alignment
            }

            // quick scan for k=0 (exact match)
            {
                // Search forward strand
                FMQuickScan scanF = FMQuickScan.scanMismatchLocations(fmIndex, q[0], Strand.FORWARD);
                if (scanF.numMismatches == 0) {
                    minMismatches = 0;
                    reportExactMatchAlignment(scanF);
                    return;
                }
                // Search reverse strand
                FMQuickScan scanR = FMQuickScan.scanMismatchLocations(fmIndex, q[1], Strand.REVERSE);
                if (scanR.numMismatches == 0) {
                    minMismatches = 0;
                    reportExactMatchAlignment(scanR);
                    return;
                }

                if (_logger.isTraceEnabled()) {
                    _logger.trace(SilkLens.toSilk("scanF", scanF));
                    _logger.trace(SilkLens.toSilk("scanR", scanR));
                }

                if (k == 0)
                    return;

                this.queryMask[0] = new QueryMask(q[0]);
                this.queryMask[1] = new QueryMask(q[1]);

                SearchState sF = null;
                SearchState sR = null;
                // Add states for both strands
                if (scanF.numMismatches <= k) {
                    if (scanF.longestMatch.start != 0 && scanF.longestMatch.start < m) {
                        // add bidirectional search state
                        sF = new SearchState(null, new Cursor(Strand.FORWARD, SearchDirection.BidirectionalForward, m,
                                scanF.longestMatch.start, scanF.longestMatch.start, null), scanF.numMismatches);
                    }
                    else {
                        sF = new SearchState(null, new Cursor(Strand.FORWARD, SearchDirection.Forward, m, 0, 0, null),
                                scanF.numMismatches);
                    }
                }

                if (scanR.numMismatches <= k) {
                    if (scanF.numMismatches > scanR.numMismatches) {
                        if (scanR.longestMatch.start != 0 && scanR.longestMatch.start < m) {
                            // add bidirectional search state
                            sR = new SearchState(null, new Cursor(Strand.REVERSE, SearchDirection.BidirectionalForward,
                                    m, scanR.longestMatch.start, scanR.longestMatch.start, null), scanR.numMismatches);
                        }
                    }
                    else {
                        sR = new SearchState(null, new Cursor(Strand.REVERSE, SearchDirection.Forward, m, 0, 0, null),
                                scanR.numMismatches);
                    }
                }

                if (scanF.numMismatches != scanR.numMismatches) {
                    if (sR != null && scanF.numMismatches < scanR.numMismatches)
                        sR.lowerThePrioity(m);
                    else if (sF != null)
                        sF.lowerThePrioity(m);
                }

                if (sF != null)
                    queue.add(sF);
                if (sR != null)
                    queue.add(sR);
            }

            while (!queue.isEmpty()) {
                SearchState c = queue.poll();
                if (_logger.isTraceEnabled()) {
                    _logger.trace("state: %s", c);
                }

                if (c.isFinished())
                    continue;

                int nm = c.getLowerBoundOfK();
                if (nm > minMismatches) {
                    ++numCutOff;
                    continue;
                }

                if (c.hasHit() || c.cursor.getRemainingBases() == 0) {
                    reportAlignment(c);
                    continue;
                }

                int allowedMismatches = k - nm;
                if (allowedMismatches < 0)
                    continue;

                if (allowedMismatches == 0) {
                    // do exact match
                    SearchState matchState = exactMatch(c);
                    if (matchState != null) {
                        reportAlignment(matchState);
                    }
                    continue;
                }

                final int strandIndex = c.cursor.getStrand().index;

                // Step forward the cursor
                // Match 
                ACGT nextBase = c.cursor.nextACGT(q);
                {
                    if (!c.isChecked(nextBase)) {
                        c.updateFlag(nextBase);
                        if (!c.siTable.isEmpty(nextBase)) {

                            SiSet nextSi = next(c, nextBase);
                            ++numFMIndexSearches;
                            SearchState nextState = c.nextState(nextBase, nextSi, queryMask[strandIndex]);
                            if (nextState != null) {
                                if (nextState.hasHit()) {
                                    reportAlignment(nextState);
                                }
                                else
                                    queue.add(nextState);
                            }
                            else
                                ++numFiltered;
                            continue; // A match is found. Proceed to the next base 
                        }
                        if (staircaseFilter.getStaircaseMask(nm + 1).get(c.cursor.getNextACGTIndex())) {
                            // mismatch is allowed at this position
                            c.lowerThePrioity(1);
                            queue.add(c); // preserve the state for back-tracking, but with lower priority
                        }
                        else {
                            numFiltered++;
                        }
                    }
                }

                // Traverse the suffix arrays for all of A, C, G and T                 
                // Add states for every bases
                for (ACGT ch : ACGT.exceptN) {
                    if (ch == nextBase)
                        continue;

                    if (!c.isChecked(ch)) {
                        c.updateFlag(ch);

                        if (!c.siTable.isEmpty(ch)) {
                            SiSet nextSi = next(c, ch);
                            ++numFMIndexSearches;
                            SearchState nextState = c.nextState(ch, nextSi, queryMask[strandIndex]);
                            if (nextState != null)
                                queue.add(nextState);
                            else
                                ++numFiltered;
                        }
                    }
                }

                // Split alignment
                c.updateSplitFlag();
                if (!c.cursor.hasSplit() && config.numSplitAlowed > 0 && nm + 1 <= minMismatches) {
                    int index = c.cursor.getIndex();
                    if (index > config.indelEndSkip && m - index > config.indelEndSkip) {
                        SearchState nextState = c.nextStateAfterSplit(nextBase);
                        if (nextState != null)
                            queue.add(nextState);
                        else
                            ++numFiltered;
                    }
                }
            }

        }

        private void reportExactMatchAlignment(FMQuickScan f) throws Exception {
            if (f.si.isUniqueHit()) {
                PosOnGenome pos = fmIndex.toGenomeCoordinate(f.si.lowerBound, m, f.strand);
                CIGAR cigar = new CIGAR();
                cigar.add(m, CIGAR.Type.Matches);
                // TODO paired-end support
                ACGTSequence r = read.getRead(0);
                String qual = read.getQual(0);
                if (!f.strand.isForward()) {
                    r = r.reverseComplement();
                    qual = reverse(qual);
                }

                AlignmentRecord rec = new AlignmentRecord(read.name(), pos.chr, f.strand, pos.pos, pos.pos + m, 0,
                        cigar, r.toString(), qual, m * config.matchScore, null);
                out.emit(rec);
            }
            else {
                // multiple hits

            }
        }

        public ReadHit verify(SearchState s, int fragmentLength, boolean isSplit) {

            SuffixInterval si = s.currentSi;
            int nm = s.getLowerBoundOfK();
            Cursor cursor = s.cursor;

            // Verification phase
            if (si == null)
                return null;

            ReadHit hit = null;
            if (si.isUniqueHit()) {
                long seqIndex = fmIndex.toCoordinate(si.lowerBound, cursor.getStrand(), cursor.getSearchDirection());
                int offset = cursor.getOffsetFromSearchHead(isSplit);
                long x = seqIndex - offset;

                if (x < 0 || x + fragmentLength > fmIndex.textSize()) {
                    return null; // ignore the match at cycle boundary
                }

                ACGTSequence ref = reference.subSequence(x, x + fragmentLength);
                ACGTSequence query = q[cursor.getStrandIndex()];
                if (isSplit) {
                    query = query.subSequence(cursor.cursorB, cursor.cursorB + fragmentLength);
                }
                if (cursor.getStrand() == Strand.REVERSE)
                    query = query.reverse();

                SWResult alignment = BitParallelSmithWaterman.alignBlock(ref, query, nm);
                hit = new ReadHit(x, fragmentLength, alignment.diff, cursor.getStrand(), null);
            }
            else {
                // TODO multiple hits

            }

            return hit;
        }

        private void reportAlignment(SearchState c) throws Exception {

            // TODO split alignment. Issue 18
            if (c.split != null) {
                // split alignment
                int splitLen = c.split.cursor.getProcessedBases();
                int m1 = Math.max(c.cursor.getProcessedBases(), m) - splitLen;
                ReadHit alignment = verify(c, m1, true);
                ReadHit splitAlignment = verify(c.split, splitLen, true);
                if (alignment == null || splitAlignment == null)
                    return;

                int diff = alignment.diff + splitAlignment.diff + 1;
                c.setLowerBoundOfK(diff);
                if (c.getLowerBoundOfK() < minMismatches)
                    minMismatches = c.getLowerBoundOfK();

                if (alignment.pos < splitLen)
                    out.emit(alignment.addSplit(splitAlignment));
                else
                    out.emit(splitAlignment.addSplit(alignment));
            }
            else {
                // Verification phase
                ReadHit alignment = verify(c, c.cursor.getFragmentLength(), false);
                if (alignment == null)
                    return; // no match

                // update the lower bound of mismatches
                c.setLowerBoundOfK(alignment.diff);

                if (c.getLowerBoundOfK() < minMismatches) {
                    minMismatches = c.getLowerBoundOfK();
                }
                out.emit(alignment);
            }

        }

        private int numFMIndexSearches  = 0;
        private int numCutOff           = 0;
        private int numFiltered         = 0;
        private int numExactSearchCount = 0;

        private SearchState exactMatch(SearchState c) {
            ++numExactSearchCount;
            Cursor cursor = c.cursor;
            final int n = cursor.getRemainingBases();
            int numExtend = 0;

            SiSet siSet = c.siTable;
            ACGT ch = c.currentACGT();
            while (numExtend < n) {
                ch = cursor.nextACGT(q);
                if (siSet.isEmpty(ch))
                    return null;
                siSet = next(c, ch);
                ++numFMIndexSearches;
                cursor = cursor.next();
                ++numExtend;
            }
            return new SearchState(null, ch, cursor, siSet, null, true, c.getLowerBoundOfK(), c.getPriority(), c.split);
        }

    }

    private static String reverse(String s) {
        if (s == null)
            return null;
        StringBuilder out = new StringBuilder(s.length());
        for (int i = s.length() - 1; i >= 0; --i)
            out.append(s.charAt(i));
        return out.toString();
    }

    /**
     * Candidate genome position of the alignment
     * 
     * @author leo
     * 
     */
    private static class CandidatePos
    {

        public final long         pos;
        public final int          k;
        public final CandidatePos chain;

        public CandidatePos(long pos, int k, CandidatePos chain) {
            this.pos = pos;
            this.k = k;
            this.chain = chain;
        }

        public static CandidatePos createFrom(SearchState s, FMIndexOnGenome fmIndex) {
            SuffixInterval si = s.currentSi;
            Cursor cursor = s.cursor;

            long seqIndex = fmIndex.toCoordinate(si.lowerBound, cursor.getStrand(), cursor.getSearchDirection());
            int offset = cursor.getCursorRange();
            long x = seqIndex - offset;

            return null;
        }

    }

    private static class Verifier
    {

    }

}
