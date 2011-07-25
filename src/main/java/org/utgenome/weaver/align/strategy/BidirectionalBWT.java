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

import java.util.Comparator;
import java.util.HashMap;
import java.util.PriorityQueue;

import org.utgenome.weaver.align.ACGT;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.AlignmentSA;
import org.utgenome.weaver.align.AlignmentScoreConfig;
import org.utgenome.weaver.align.BitVector;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.Range;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.SuffixInterval;
import org.utgenome.weaver.align.record.AlignmentRecord;
import org.utgenome.weaver.align.record.RawRead;
import org.utgenome.weaver.align.record.ReadSequence;
import org.utgenome.weaver.align.strategy.BidirectionalBWT.Alignment.ExtensionType;
import org.utgenome.weaver.parallel.Reporter;
import org.xerial.util.ObjectHandlerBase;

/**
 * Alignment algorithm using Bi-directional BWT
 * 
 * @author leo
 * 
 */
public class BidirectionalBWT
{
    private final FMIndexOnGenome      fmIndex;
    private final Reporter             reporter;
    private final long                 N;
    private final AlignmentScoreConfig config;

    public BidirectionalBWT(FMIndexOnGenome fmIndex, Reporter reporter) {
        this.fmIndex = fmIndex;
        this.reporter = reporter;
        this.N = fmIndex.textSize();
        this.config = new AlignmentScoreConfig();
    }

    public static class SARange
    {
        public final SuffixInterval si;
        public final Strand         strand;

        public SARange(SuffixInterval si, Strand strand) {
            this.si = si;
            this.strand = strand;
        }
    }

    public static class QuickScanResult
    {
        public final SuffixInterval si;
        public final BitVector      breakPoint;
        public final int            numMismatches;
        public final Range          longestMatch;

        public QuickScanResult(SuffixInterval si, BitVector breakPoint, int numMismatches, Range longestMatch) {
            this.si = si;
            this.breakPoint = breakPoint;
            this.numMismatches = numMismatches;
            this.longestMatch = longestMatch;
        }
    }

    public QuickScanResult scanMismatchLocations(ACGTSequence query, Strand strand) {
        int qLen = (int) query.textSize();
        int numMismatches = 0;
        BitVector breakPoint = new BitVector(qLen);
        SuffixInterval si = new SuffixInterval(0, N - 1);
        int longestMatchLength = 0;
        int mark = 0;
        Range longestMatch = null;
        int i = 0;
        for (; i < qLen; ++i) {
            ACGT ch = query.getACGT(i);
            si = fmIndex.forwardSearch(strand, ch, si);
            if (!si.isValidRange()) {
                si = new SuffixInterval(0, N - 1);
                breakPoint.set(i, true);
                numMismatches++;
                if (longestMatch == null || longestMatch.length() < (i - mark)) {
                    longestMatch = new Range(mark, i);
                }
                mark = i + 1;
            }
        }
        if (longestMatch == null || longestMatch.length() < (i - mark)) {
            longestMatch = new Range(mark, i);
        }

        return new QuickScanResult(si, breakPoint, numMismatches, longestMatch);
    }

    void report(AlignmentSA result) throws Exception {
        fmIndex.toGenomeCoordinate(result, new ObjectHandlerBase<AlignmentRecord>() {
            @Override
            public void handle(AlignmentRecord r) throws Exception {
                reporter.emit(r);
            }
        });
    }

    void report(Alignment result) {

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

    public Alignment prepareInitialAlignmentState(ACGTSequence q, QuickScanResult scan, Strand strand) {
        // Break the read sequence into three parts 
        int s1 = (int) q.textSize() / 3;
        int s2 = (int) q.textSize() / 3 * 2;

        int[] m = new int[3];
        m[0] = scan.breakPoint.countOneBits(0L, s1);
        m[1] = scan.breakPoint.countOneBits(s1, s2);
        m[2] = scan.breakPoint.countOneBits(s2, q.textSize());

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
            return ForwardAlignment.newInstance(q, strand, N);
        case MIDDLE:
            return BidirectionalForwardAlignment.newInstance(q, strand, N);
        case TAIL:
            return BackwardAlignment.newInstance(q, strand, N);
        }
        return null;
    }

    public static class AlignmentQueue
    {
        private final PriorityQueue<Alignment> queue;
        private final AlignmentScoreConfig     config;
        private int                            scoreLowerBound;

        public AlignmentQueue(AlignmentScoreConfig config) {
            this.config = config;
            this.queue = new PriorityQueue<Alignment>(11, new Comparator<Alignment>() {
                @Override
                public int compare(Alignment o1, Alignment o2) {
                    // If the upper bound of the score is larger than the other, search it first
                    int diff = o2.getUpperBoundOfScore(AlignmentQueue.this.config)
                            - o1.getUpperBoundOfScore(AlignmentQueue.this.config);
                    if (diff != 0)
                        return diff;

                    return o2.getRemainingBases() - o1.getRemainingBases();
                }
            });
        }

        @Override
        public String toString() {
            return String.format("size:%d, lb:%d", queue.size(), scoreLowerBound);
        }

        public Alignment poll() {
            return queue.poll();
        }

        public boolean isEmpty() {
            return queue.isEmpty();
        }

        public boolean add(Alignment e) {
            if (e.score > scoreLowerBound)
                scoreLowerBound = e.score;

            int upperBound = e.getUpperBoundOfScore(config);

            if (upperBound < scoreLowerBound) {
                // Discard the alignment whose score cannot exceed the lower bound
                return false;
            }

            return queue.add(e);
        }

    }

    public void align(RawRead r) throws Exception {

        // TODO PE mapping
        ReadSequence read = (ReadSequence) r;

        ACGTSequence qF = new ACGTSequence(read.seq);

        // Find potential mismatch positions for forward direction
        QuickScanResult scanF = scanMismatchLocations(qF, Strand.FORWARD);
        if (scanF.numMismatches == 0) {
            // Found an exact match
            report(AlignmentSA.exactMatch(config, r.name(), qF, scanF.si, Strand.FORWARD));
            return;
        }

        // Find potential mismatch positions for reverse direction
        ACGTSequence qC = qF.complement();
        QuickScanResult scanR = scanMismatchLocations(qC, Strand.REVERSE);
        if (scanR.numMismatches == 0) {
            // Found an exact match
            report(AlignmentSA.exactMatch(config, r.name(), qC, scanR.si, Strand.REVERSE));
            return;
        }

        AlignmentQueue alignmentQueue = new AlignmentQueue(config);
        // Set the initial search states
        alignmentQueue.add(prepareInitialAlignmentState(qF, scanF, Strand.FORWARD));
        alignmentQueue.add(prepareInitialAlignmentState(qC, scanR, Strand.REVERSE));

        AlignmentVisitor visitor = new AlignmentProcessor(fmIndex, alignmentQueue, config);

        // Search iteration
        while (!alignmentQueue.isEmpty()) {
            Alignment current = alignmentQueue.poll();

            if (current.isFinished()) {
                report(current);
            }

            // extend the match
            current.accept(visitor);
        }

    }

    public static abstract class Alignment
    {
        public static enum ExtensionType {
            MATCH, INSERTION, DELETION
        };

        public final ACGTSequence  read;
        public final Strand        strand;
        public final ExtensionType extensionType;
        public final int           cursor;
        public final int           score;
        public final int           numMismatches;
        public final int           numGapOpens;
        public final int           numGapExtend;

        protected Alignment(ACGTSequence read, Strand strand, ExtensionType extensionType, int cursor, int score,
                int numMismatches, int numGapOpens, int numGapExtend) {
            this.read = read;
            this.strand = strand;
            this.extensionType = extensionType;
            this.cursor = cursor;
            this.score = score;
            this.numMismatches = numMismatches;
            this.numGapOpens = numGapOpens;
            this.numGapExtend = numGapExtend;
        }

        protected Alignment(Alignment other) {
            this.read = other.read;
            this.strand = other.strand;
            this.extensionType = other.extensionType;
            this.cursor = other.cursor;
            this.score = other.score;
            this.numMismatches = other.numMismatches;
            this.numGapOpens = other.numGapOpens;
            this.numGapExtend = other.numGapExtend;
        }

        @Override
        public String toString() {
            return String.format("strand:%s, score:%d, cursor:%d, state:%s, mm:%d, go:%d, ge:%d", strand.symbol, score,
                    cursor, extensionType, numMismatches, numGapOpens, numGapExtend);
        }

        public abstract int getRemainingBases();

        public abstract int getUpperBoundOfScore(AlignmentScoreConfig config);

        public abstract boolean isFinished();

        public abstract void accept(AlignmentVisitor visitor);

    }

    public static interface AlignmentVisitor
    {
        public void forwardAlignment(ForwardAlignment forward);

        public void bidirectionalForwardAlignment(BidirectionalForwardAlignment biForward);

        public void backwardAlignment(BackwardAlignment forward);
    }

    public static class AlignmentProcessor implements AlignmentVisitor
    {
        private final FMIndexOnGenome      fmIndex;
        private final AlignmentQueue       queue;
        private final AlignmentScoreConfig config;

        public AlignmentProcessor(FMIndexOnGenome fmIndex, AlignmentQueue queue, AlignmentScoreConfig config) {
            this.fmIndex = fmIndex;
            this.queue = queue;
            this.config = config;
        }

        public void exactMatchForward(ForwardAlignment c) {
            int cursor = c.cursor;
            int score = c.score;
            SuffixInterval si = c.reverseSi;
            while (cursor < c.read.textSize()) {
                si = fmIndex.forwardSearch(c.strand, c.read.getACGT(cursor), si);
                if (!si.isValidRange())
                    return;
                cursor++;
                score += config.matchScore;
            }
            queue.add(new ForwardAlignment(c.read, c.strand, ExtensionType.MATCH, cursor, score, c.numMismatches,
                    c.numGapOpens, c.numGapExtend, si));
        }

        @Override
        public void forwardAlignment(ForwardAlignment c) {
            int cursor = c.cursor;

            SuffixInterval si = c.reverseSi;
            if (cursor >= c.read.textSize()) {
                return;
            }

            if (c.getUpperBoundOfScore(config) < queue.scoreLowerBound)
                return; // no need to proceed

            if (c.numMismatches >= config.numMismatchesAllowed) {
                // exact match
                exactMatchForward(c);
                return;
            }

            // Compute the next SA ranges for A, C, G, T, N
            SuffixInterval[] next = new SuffixInterval[ACGT.values().length]; // for A, C, G, T, N
            {
                int i = 0;
                for (ACGT ch : ACGT.values()) {
                    next[i++] = fmIndex.forwardSearch(c.strand, ch, si);
                }
            }

            // Search for indels
            switch (c.extensionType) {
            case MATCH: { // gap open
                if (c.numGapOpens < config.numGapOpenAllowed) {
                    // insertion to reference
                    queue.add(new ForwardAlignment(c.read, c.strand, ExtensionType.INSERTION, c.cursor + 1, c.score
                            - config.gapOpenPenalty, c.numMismatches, c.numGapOpens + 1, c.numGapExtend, c.reverseSi));
                    // deletion from reference
                    for (ACGT ch : ACGT.exceptN) {
                        if (next[ch.code].isValidRange())
                            queue.add(new ForwardAlignment(c.read, c.strand, ExtensionType.DELETION, c.cursor, c.score
                                    - config.gapOpenPenalty, c.numMismatches, c.numGapOpens + 1, c.numGapExtend,
                                    next[ch.code]));
                    }
                    break;
                }
            }
            case DELETION: { // gap extension
                for (ACGT ch : ACGT.exceptN) {
                    if (next[ch.code].isValidRange())
                        queue.add(new ForwardAlignment(c.read, c.strand, ExtensionType.DELETION, c.cursor, c.score
                                - config.gapExtentionPenalty, c.numMismatches, c.numGapOpens, c.numGapExtend + 1,
                                next[ch.code]));
                }
                break;
            }
            case INSERTION: { // gap extension
                queue.add(new ForwardAlignment(c.read, c.strand, ExtensionType.INSERTION, c.cursor + 1, c.score
                        - config.gapExtentionPenalty, c.numMismatches, c.numGapOpens, c.numGapExtend + 1, c.reverseSi));
                break;
            }
            }

            // Search for mismatches
            if (c.numMismatches < config.numMismatchesAllowed) {
                for (ACGT ch : ACGT.exceptN) {
                    SuffixInterval nextSi = next[ch.code];
                    if (!nextSi.isValidRange())
                        continue;

                    if (ch == c.read.getACGT(c.cursor)) {
                        // match
                        queue.add(new ForwardAlignment(c.read, c.strand, ExtensionType.MATCH, c.cursor + 1, c.score
                                + config.matchScore, c.numMismatches, c.numGapOpens, c.numGapExtend, nextSi));
                    }
                    else {
                        // mismatch
                        queue.add(new ForwardAlignment(c.read, c.strand, ExtensionType.MATCH, c.cursor + 1, c.score
                                - config.mismatchPenalty, c.numMismatches + 1, c.numGapOpens, c.numGapExtend, nextSi));
                    }
                }
            }
            else {
                // exact match only
                queue.add(new ForwardAlignment(c.read, c.strand, ExtensionType.MATCH, c.cursor + 1, c.score
                        + config.matchScore, c.numMismatches, c.numGapOpens, c.numGapExtend, next[c.read
                        .getACGT(c.cursor).code]));
            }

        }

        @Override
        public void bidirectionalForwardAlignment(BidirectionalForwardAlignment biForward) {
            // TODO Auto-generated method stub

        }

        @Override
        public void backwardAlignment(BackwardAlignment forward) {
            // TODO Auto-generated method stub

        }

    }

    public static class ForwardAlignment extends Alignment
    {
        public final SuffixInterval reverseSi;

        public static ForwardAlignment newInstance(ACGTSequence read, Strand strand, long N) {
            return new ForwardAlignment(read, strand, ExtensionType.MATCH, 0, 0, 0, 0, 0, new SuffixInterval(0, N));
        }

        public ForwardAlignment(ACGTSequence read, Strand strand, ExtensionType extensionType, int cursor, int score,
                int numMismatches, int numGapOpens, int numGapExtend, SuffixInterval reverseSi) {
            super(read, strand, extensionType, cursor, score, numMismatches, numGapOpens, numGapExtend);
            this.reverseSi = reverseSi;
        }

        public ForwardAlignment(Alignment base, SuffixInterval reverseSi) {
            super(base);
            this.reverseSi = reverseSi;
        }

        @Override
        public String toString() {
            return String.format("%s%d:%d %d/%d/%d/%s", extensionType.name().charAt(0), cursor, score, numMismatches,
                    numGapOpens, numGapExtend, strand.symbol);
        }

        @Override
        public int getRemainingBases() {
            return Math.max(0, (int) read.textSize() - cursor);
        }

        @Override
        public int getUpperBoundOfScore(AlignmentScoreConfig config) {
            return score + config.matchScore * getRemainingBases();
        }

        @Override
        public boolean isFinished() {
            return cursor >= read.textSize();
        }

        @Override
        public void accept(AlignmentVisitor visitor) {
            visitor.forwardAlignment(this);
        }

    }

    public static class BackwardAlignment extends Alignment
    {
        public final SuffixInterval forwardSi;

        public static BackwardAlignment newInstance(ACGTSequence read, Strand strand, long N) {
            return new BackwardAlignment(read, strand, ExtensionType.MATCH, (int) read.textSize() - 1, 0, 0, 0, 0,
                    new SuffixInterval(0, N));
        }

        public BackwardAlignment(ACGTSequence read, Strand strand, ExtensionType extensionType, int cursor, int score,
                int numMismatches, int numGapOpens, int numGapExtend, SuffixInterval forwardSi) {
            super(read, strand, extensionType, cursor, score, numMismatches, numGapOpens, numGapExtend);
            this.forwardSi = forwardSi;
        }

        @Override
        public int getRemainingBases() {
            return Math.max(0, cursor + 1);
        }

        @Override
        public int getUpperBoundOfScore(AlignmentScoreConfig config) {
            return score + config.matchScore * getRemainingBases();
        }

        @Override
        public boolean isFinished() {
            return cursor < 0;
        }

        @Override
        public void accept(AlignmentVisitor visitor) {
            visitor.backwardAlignment(this);
        }

    }

    public static class BidirectionalForwardAlignment extends Alignment
    {
        public final SuffixInterval forwardSi;
        public final SuffixInterval backwardSi;
        public final int            reverseCursor;

        public static BidirectionalForwardAlignment newInstance(ACGTSequence read, Strand strand, long N) {
            final int K = (int) read.textSize();
            int s2 = K / 3;
            return new BidirectionalForwardAlignment(read, strand, ExtensionType.MATCH, s2, 0, 0, 0, 0,
                    new SuffixInterval(0, N), new SuffixInterval(0, N), s2);
        }

        protected BidirectionalForwardAlignment(ACGTSequence read, Strand strand, ExtensionType extensionType,
                int cursor, int score, int numMismatches, int numGapOpens, int numGapExtend, SuffixInterval forwardSi,
                SuffixInterval backwardSi, int reverseCursor) {
            super(read, strand, extensionType, cursor, score, numMismatches, numGapOpens, numGapExtend);
            this.forwardSi = forwardSi;
            this.backwardSi = backwardSi;
            this.reverseCursor = reverseCursor;
        }

        @Override
        public int getRemainingBases() {
            int remainingLeft = (int) read.textSize() / 3;
            int remainingRight = (int) read.textSize() - cursor;
            int remainingBases = remainingLeft + remainingRight;
            return remainingBases;
        }

        @Override
        public int getUpperBoundOfScore(AlignmentScoreConfig config) {
            return score + config.matchScore * getRemainingBases();
        }

        @Override
        public boolean isFinished() {
            return false;
        }

        @Override
        public void accept(AlignmentVisitor visitor) {
            visitor.bidirectionalForwardAlignment(this);
        }
    }

}
