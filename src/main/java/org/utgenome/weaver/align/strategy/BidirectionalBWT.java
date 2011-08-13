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
import org.utgenome.weaver.align.AlignmentSA;
import org.utgenome.weaver.align.AlignmentScoreConfig;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.SuffixInterval;
import org.utgenome.weaver.align.record.AlignmentRecord;
import org.utgenome.weaver.align.record.RawRead;
import org.utgenome.weaver.parallel.Reporter;
import org.xerial.lens.SilkLens;
import org.xerial.util.ObjectHandlerBase;
import org.xerial.util.log.Logger;

/**
 * Alignment algorithm using Bi-directional BWT.
 * 
 * @see SuffixFilter for better integration with bidirectional BWT and staircase
 *      NFA
 * 
 * @author leo
 * 
 */
public class BidirectionalBWT
{
    private static Logger         _logger = Logger.getLogger(BidirectionalBWT.class);

    private final FMIndexOnGenome fmIndex;
    private final Reporter        reporter;
    private final long            N;
    private AlignmentScoreConfig  config;

    public BidirectionalBWT(FMIndexOnGenome fmIndex, Reporter reporter) {
        this.fmIndex = fmIndex;
        this.reporter = reporter;
        this.N = fmIndex.textSize();
        this.config = new AlignmentScoreConfig();
    }

    public void setAlignmentScoreConfig(AlignmentScoreConfig config) {
        this.config = config;
    }

    void report(AlignmentSA result) throws Exception {
        fmIndex.toGenomeCoordinate(result, new ObjectHandlerBase<AlignmentRecord>() {
            @Override
            public void handle(AlignmentRecord r) throws Exception {
                reporter.emit(r);
            }
        });
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

    public BWAState prepareInitialAlignmentState(ACGTSequence q, QuickScanResult scan, Strand strand) {
        // Break the read sequence into three parts
        int M = (int) q.textSize();
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
        case FRONT:
            return new BWAState(q, strand, SearchDirection.Forward, ExtensionType.MATCH, 0, -1, Score.initial(),
                    fmIndex.wholeSARange());
        case MIDDLE: {
            return new BWAState(q, strand, SearchDirection.BidirectionalForward, ExtensionType.MATCH, s1, s1,
                    Score.initial(), fmIndex.wholeSARange());
        }
        case TAIL:
            return new BWAState(q, strand, SearchDirection.Backward, ExtensionType.MATCH, M, M, Score.initial(),
                    fmIndex.wholeSARange());
        }
        return null;
    }

    public BWAState exactMatch(BWAState aln) {
        while (!aln.isFinished()) {
            SuffixInterval nextSi = aln.nextSi(fmIndex, aln.nextACGT(), aln.si);
            if (nextSi.isEmpty())
                return null;
            aln = aln.extendWithMatch(config, nextSi);
        }
        return aln;
    }

    public void align(RawRead r) throws Exception {

        // TODO PE mapping
        ACGTSequence qF = r.getRead(0);
        _logger.debug("query: " + qF);

        if (qF.fastCount(ACGT.N, 0, qF.textSize()) > config.maximumEditDistances) {
            // too many Ns in the query sequence
            return;
        }

        // Find potential mismatch positions for forward direction
        QuickScanResult scanF = QuickScanResult.scanMismatchLocations(fmIndex, qF, Strand.FORWARD);
        if (scanF.numMismatches == 0) {
            // Found an exact match
            report(AlignmentSA.exactMatch(config, r.name(), qF, scanF.si, Strand.FORWARD));
            return;
        }

        // Find potential mismatch positions for reverse direction
        ACGTSequence qC = qF.complement();
        QuickScanResult scanR = QuickScanResult.scanMismatchLocations(fmIndex, qC, Strand.REVERSE);
        if (scanR.numMismatches == 0) {
            // Found an exact match
            report(AlignmentSA.exactMatch(config, r.name(), qC, scanR.si, Strand.REVERSE));
            return;
        }

        if (_logger.isDebugEnabled()) {
            _logger.debug(SilkLens.toSilk("scanF", scanF));
            _logger.debug(SilkLens.toSilk("scanR", scanR));
        }

        AlignmentQueue queue = new AlignmentQueue(config);
        // Set the initial search states
        queue.add(prepareInitialAlignmentState(qF, scanF, Strand.FORWARD));
        queue.add(prepareInitialAlignmentState(qC, scanR, Strand.REVERSE));

        // Search iteration
        while (!queue.isEmpty()) {
            BWAState c = queue.poll(); // current 

            if (c.isFinished() && c.score.score >= queue.bestScore) {
                report(c);
                queue.bestScore = c.score.score;
                continue;
            }

            SuffixInterval si = c.suffixInterval();
            if (c.isFinished()) {
                continue;
            }

            int upperBound = c.getUpperBoundOfScore(config);

            if (upperBound < queue.bestScore)
                continue; // no need to proceed

            int remainingDist = config.maximumEditDistances
                    - (c.score.numMismatches + c.score.numGapOpens + c.score.numGapExtend);
            if (remainingDist < 0)
                continue;

            if (_logger.isDebugEnabled())
                _logger.debug("[%3d] %s SI:%s ", queue.queue.size(), c, c.suffixInterval());

            if (remainingDist == 0) {
                if (c.extensionType == ExtensionType.MATCH) {
                    // exact match
                    BWAState a = exactMatch(c);
                    if (a != null)
                        queue.add(a);
                }
                continue;
            }

            // Compute the next SA ranges for A, C, G, T, N
            SuffixInterval[] next = new SuffixInterval[ACGT.values().length]; // for A, C, G, T, N
            {
                int i = 0;
                for (ACGT ch : ACGT.values()) {
                    next[i++] = c.nextSi(fmIndex, ch, si);
                }
            }

            // Search for indels
            switch (c.extensionType) {
            case MATCH: { // gap open
                if (c.gapOpenIsAllowed(config) && upperBound - config.gapOpenPenalty > queue.bestScore) {
                    // insertion to reference
                    queue.add(c.startInsertion(config));
                    // deletion from reference
                    for (ACGT ch : ACGT.exceptN) {
                        if (!next[ch.code].isEmpty())
                            queue.add(c.startDeletion(config, next[ch.code]));
                    }
                    break;
                }
            }
            case DELETION: { // gap extension
                if (c.gapExtensionIsAllowed(config) && upperBound - config.gapExtensionPenalty > queue.bestScore) {
                    for (ACGT ch : ACGT.exceptN) {
                        if (!next[ch.code].isEmpty())
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
            for (ACGT ch : ACGT.exceptN) {
                SuffixInterval nextSi = next[ch.code];
                if (nextSi.isEmpty())
                    continue;

                if (ch == c.nextACGT()) {
                    // match
                    queue.add(c.extendWithMatch(config, nextSi));
                }
                else {
                    if (remainingDist > 0) {
                        // mismatch
                        queue.add(c.extendWithMisMatch(config, nextSi));
                    }
                }
            }

        }
        if (_logger.isDebugEnabled())
            _logger.debug("push count: %,d", queue.pushCount);

    }

}
