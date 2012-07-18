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
// BidirectionalBWT.java
// Since: 2011/06/29
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.strategy;

import java.util.HashMap;

import org.utgenome.weaver.align.ACGT;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.Aligner;
import org.utgenome.weaver.align.AlignmentConfig;
import org.utgenome.weaver.align.AlignmentScoreConfig;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.SiSet;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.record.Read;
import org.utgenome.weaver.parallel.Reporter;
import org.xerial.lens.SilkLens;
import org.xerial.util.log.Logger;

/**
 * Alignment algorithm using Bi-directional BWT.
 * 
 * @see BidirectionalSuffixFilter for better integration with bidirectional BWT
 *      and staircase NFA
 * 
 * @author leo
 * 
 */
public class BidirectionalBWT implements Aligner
{
    private static Logger         _logger                     = Logger.getLogger(BidirectionalBWT.class);

    private final FMIndexOnGenome fmIndex;
    private Reporter              reporter;
    private final long            N;
    private AlignmentConfig       config;
    private boolean               disableBidirecdtionalSearch = false;

    public BidirectionalBWT(FMIndexOnGenome fmIndex, Reporter reporter, AlignmentConfig config) {
        this.fmIndex = fmIndex;
        this.reporter = reporter;
        this.N = fmIndex.textSize();
        this.config = config;
    }

    public void disableBidirectionalSearch() {
        this.disableBidirecdtionalSearch = true;
    }

    void report(BWAState result) throws Exception {
        reporter.emit(result);
    }

    private static enum SearchStart {
        FRONT, MIDDLE, TAIL
    }

    private static HashMap<Integer, SearchStart> searchStrategyTable = new HashMap<Integer, SearchStart>();

    static {
        searchStrategyTable.put(Integer.parseInt("001", 2), SearchStart.FRONT);
        searchStrategyTable.put(Integer.parseInt("010", 2), SearchStart.FRONT);
        searchStrategyTable.put(Integer.parseInt("100", 2), SearchStart.TAIL);
        searchStrategyTable.put(Integer.parseInt("101", 2), SearchStart.MIDDLE);
        searchStrategyTable.put(Integer.parseInt("110", 2), SearchStart.TAIL);
        searchStrategyTable.put(Integer.parseInt("011", 2), SearchStart.FRONT);
    }

    public BWAState prepareInitialAlignmentState(ACGTSequence q, FMQuickScan scan, Strand strand) {

        int M = (int) q.textSize();

        if (disableBidirecdtionalSearch) {
            return new BWAState(new Cursor(strand, SearchDirection.Forward, 0, M, 0, 0), ExtensionType.MATCH,
                    Score.initial(), fmIndex.initSet(SearchDirection.Forward));
        }

        // Break the read sequence into three parts

        int s1 = M / 3;
        int s2 = M / 3 * 2;

        int[] m = new int[3];
        m[0] = scan.breakPoint.countOneBits(0L, s1);
        m[1] = scan.breakPoint.countOneBits(s1, s2);
        m[2] = scan.breakPoint.countOneBits(s2, M);

        int mismatchBlockFlag = 0;
        int flag = 4;
        for (int i = 0; i < 3; ++i) {
            if (m[i] > 0)
                mismatchBlockFlag += flag;
            flag >>>= 1;
        }

        // Get the search start location according to the mismatch positions
        SearchStart searchStart = searchStrategyTable.get(mismatchBlockFlag);
        if (searchStart == null) {
            // Start the search from the smallest mismatch block
            int minBlock = 0;
            for (int i = 1; i < 3; ++i) {
                if (m[i - 1] > m[i])
                    minBlock = i;
            }
            switch (minBlock) {
            case 0:
                searchStart = SearchStart.FRONT;
                break;
            case 1:
                searchStart = SearchStart.MIDDLE;
                break;
            case 2:
                searchStart = SearchStart.TAIL;
                break;
            }
        }

        switch (searchStart) {
        default:
        case FRONT:
            return new BWAState(new Cursor(strand, SearchDirection.Forward, 0, M, 0, 0), ExtensionType.MATCH,
                    Score.initial(), fmIndex.initSet(SearchDirection.Forward));
        case MIDDLE: {
            return new BWAState(new Cursor(strand, SearchDirection.BidirectionalForward, 0, M, s1, s1),
                    ExtensionType.MATCH, Score.initial(), fmIndex.initSet(SearchDirection.BidirectionalForward));
        }
        case TAIL:
            return new BWAState(new Cursor(strand, SearchDirection.Backward, 0, M, M, M), ExtensionType.MATCH,
                    Score.initial(), fmIndex.initSet(SearchDirection.Backward));
        }
    }

    public BWAState exactMatch(BWAState aln, ACGTSequence[] q) {

        Cursor cursor = aln.cursor;
        final int n = cursor.getRemainingBases();
        int numExtend = 0;

        SiSet siSet = aln.si;
        ACGT ch = null;
        while (numExtend < n) {
            ch = cursor.nextACGT(q);
            if (siSet.isEmpty(ch))
                return null;
            siSet = cursor.nextSi(fmIndex, siSet, ch);
            cursor = cursor.next();
            ++numExtend;
        }

        return new BWAState(cursor, ExtensionType.MATCH, aln.score.extendWithMatch(config, n), siSet);
    }

    public void align(Read r) throws Exception {

        // TODO PE mapping
        ACGTSequence qF = r.getRead(0);
        if (_logger.isTraceEnabled())
            _logger.trace("query: " + qF);

        if (qF.fastCount(ACGT.N, 0, qF.textSize()) > config.getMaximumEditDistance(qF.length())) {
            // too many Ns in the query sequence
            return;
        }

        // Find potential mismatch positions for forward direction
        FMQuickScan scanF = FMQuickScan.scanMismatchLocations(fmIndex, qF, Strand.FORWARD);
        if (scanF.numMismatches == 0) {
            // Found an exact match
            reporter.emit(scanF);
            return;
        }

        // Find potential mismatch positions for reverse direction
        ACGTSequence qC = qF.complement();
        FMQuickScan scanR = FMQuickScan.scanMismatchLocations(fmIndex, qC, Strand.REVERSE);
        if (scanR.numMismatches == 0) {
            // Found an exact match
            reporter.emit(scanR);
            return;
        }

        ACGTSequence[] q = new ACGTSequence[] { qF, qC };

        if (_logger.isTraceEnabled()) {
            _logger.trace(SilkLens.toSilk("scanF", scanF));
            _logger.trace(SilkLens.toSilk("scanR", scanR));
        }

        AlignmentQueue queue = new AlignmentQueue(config);
        // Set the initial search states
        queue.add(prepareInitialAlignmentState(qF, scanF, Strand.FORWARD));
        queue.add(prepareInitialAlignmentState(qC, scanR, Strand.REVERSE));

        int numFMIndexSearches = 0;
        final int m = (int) qF.textSize();
        // Search iteration
        while (!queue.isEmpty()) {
            BWAState c = queue.poll(); // current 

            if (c.isFinished() && c.score.score >= queue.bestScore) {
                report(c);
                queue.bestScore = c.score.score;
                continue;
            }

            int upperBound = c.getUpperBoundOfScore(config);

            if (upperBound < queue.bestScore)
                continue; // no need to proceed

            int posInRead = c.cursor.getProcessedBases();

            int remainingDist = config.getMaximumEditDistance(m)
                    - (c.score.numMismatches + c.score.numGapOpens + c.score.numGapExtend);
            if (remainingDist < 0)
                continue;

            if (remainingDist == 0) {
                if (c.extensionType == ExtensionType.MATCH) {
                    // exact match
                    BWAState a = exactMatch(c, q);
                    if (a != null)
                        queue.add(a);
                }
                continue;
            }

            // Compute next suffix intervals for A, C, G, T
            SiSet[] next = new SiSet[ACGT.exceptN.length];
            for (ACGT ch : ACGT.exceptN) {
                next[ch.code] = c.nextSi(fmIndex, ch);
                ++numFMIndexSearches;
            }

            // Search for indels
            switch (c.extensionType) {
            case MATCH: { // gap open
                if (c.gapOpenIsAllowed(config, posInRead, m) && upperBound - config.gapOpenPenalty > queue.bestScore) {
                    // insertion to reference
                    queue.add(c.startInsertion(config));
                    // deletion from reference
                    for (ACGT ch : ACGT.exceptN) {
                        if (!c.si.isEmpty(ch)) {
                            queue.add(c.startDeletion(config, next[ch.code]));
                        }
                    }
                    break;
                }
            }
            case DELETION: { // gap extension
                if (c.gapExtensionIsAllowed(config) && upperBound - config.gapExtensionPenalty > queue.bestScore) {
                    for (ACGT ch : ACGT.exceptN) {
                        if (!c.si.isEmpty(ch))
                            queue.add(c.extendDeletion(config, next[ch.code]));
                    }
                }
                break;
            }
            case INSERTION: { // gap extension
                if (c.gapExtensionIsAllowed(config) && upperBound - config.gapExtensionPenalty > queue.bestScore) {
                    queue.add(c.extendInsertion(config));
                }
                break;
            }
            }

            // Search for mismatches
            if (remainingDist > 0) {
                ACGT nextCh = c.nextACGT(q);
                for (ACGT ch : ACGT.exceptN) {
                    if (c.si.isEmpty(ch))
                        continue;

                    if (ch == nextCh) {
                        // match
                        queue.add(c.extendWithMatch(config, next[ch.code]));
                    }
                    else {
                        // mismatch
                        queue.add(c.extendWithMisMatch(config, next[ch.code]));
                    }
                }
            }

        }
        if (_logger.isDebugEnabled())
            _logger.debug("FM Search:%d, push count: %,d", numFMIndexSearches, queue.pushCount);

    }

    @Override
    public void align(Read read, Reporter out) throws Exception {
        this.reporter = out;
        align(read);
    }

}
