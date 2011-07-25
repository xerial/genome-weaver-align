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
        private int                            bestScore;

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

                    return o1.getRemainingBases() - o2.getRemainingBases();
                }
            });
        }

        @Override
        public String toString() {
            return String.format("size:%d, best:%d", queue.size(), bestScore);
        }

        public Alignment poll() {
            return queue.poll();
        }

        public boolean isEmpty() {
            return queue.isEmpty();
        }

        public boolean add(Alignment e) {
            if (e.getUpperBoundOfScore(config) < bestScore)
                return false;

            if (e.isFinished() && e.score.score > bestScore)
                bestScore = e.score.score;

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

            if (current.isFinished() && current.score.score >= alignmentQueue.bestScore) {
                report(current);
                alignmentQueue.bestScore = current.score.score;
                continue;
            }

            // extend the match
            current.accept(visitor);
        }

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

        public void exactMatch(Alignment c) {
            int cursor = c.cursor;
            int score = c.score.score;
            SuffixInterval si = c.suffixInterval();
            while (cursor < c.read.textSize()) {
                si = c.nextSi(fmIndex, c.getACGT(cursor), si);
                if (!si.isValidRange())
                    return;
                cursor++;
                score += config.matchScore;
            }
            queue.add(c.newInstance(score, cursor, si));
        }

        public void approximateMatch(Alignment c) {
            int cursor = c.cursor;

            SuffixInterval si = c.suffixInterval();
            if (c.isFinished()) {
                return;
            }

            if (c.getUpperBoundOfScore(config) < queue.bestScore)
                return; // no need to proceed

            if (!c.mismatchIsAllowed(config)) {
                // exact match
                exactMatch(c);
                return;
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
                if (c.gapIsAllowed(config)) {
                    // insertion to reference
                    queue.add(c.startInsertion(config));
                    // deletion from reference
                    for (ACGT ch : ACGT.exceptN) {
                        if (next[ch.code].isValidRange())
                            queue.add(c.startDeletion(config, next[ch.code]));
                    }
                    break;
                }
            }
            case DELETION: { // gap extension
                for (ACGT ch : ACGT.exceptN) {
                    if (next[ch.code].isValidRange())
                        queue.add(c.extendDeletion(config, next[ch.code]));
                }
                break;
            }
            case INSERTION: { // gap extension
                queue.add(c.extendInsertion(config));
                break;
            }
            }

            // Search for mismatches
            for (ACGT ch : ACGT.exceptN) {
                SuffixInterval nextSi = next[ch.code];
                if (!nextSi.isValidRange())
                    continue;

                if (ch == c.nextACGT()) {
                    // match
                    queue.add(c.extendWithMatch(config, nextSi));
                }
                else {
                    // mismatch
                    queue.add(c.extendWithMisMatch(config, nextSi));
                }
            }

        }

        @Override
        public void forwardAlignment(ForwardAlignment c) {
            approximateMatch(c);
        }

        @Override
        public void backwardAlignment(BackwardAlignment c) {
            approximateMatch(c);
        }

        @Override
        public void bidirectionalForwardAlignment(BidirectionalForwardAlignment biForward) {
            // TODO Auto-generated method stub

        }

    }

    public static abstract class Alignment
    {
        public static class Score
        {
            public final int score;
            public final int numMismatches;
            public final int numGapOpens;
            public final int numGapExtend;

            public Score(int score, int numMismatches, int numGapOpens, int numGapExtend) {
                this.score = score;
                this.numMismatches = numMismatches;
                this.numGapOpens = numGapOpens;
                this.numGapExtend = numGapExtend;
            }

            @Override
            public String toString() {
                return String.format("%3d:%d/%d/%d", score, numMismatches, numGapOpens, numGapExtend);
            }

            public static Score initial() {
                return new Score(0, 0, 0, 0);
            }

            public Score update(int newScore) {
                return new Score(newScore, numMismatches, numGapOpens, numGapExtend);
            }

            public Score extendWithMatch(AlignmentScoreConfig config) {
                return new Score(score + config.matchScore, numMismatches, numGapOpens, numGapExtend);
            }

            public Score extendWithMismatch(AlignmentScoreConfig config) {
                return new Score(score - config.mismatchPenalty, numMismatches + 1, numGapOpens, numGapExtend);
            }

            public Score extendWithGapOpen(AlignmentScoreConfig config) {
                return new Score(score - config.gapOpenPenalty, numMismatches, numGapOpens + 1, numGapExtend);
            }

            public Score extendWithGapExtend(AlignmentScoreConfig config) {
                return new Score(score - config.gapExtentionPenalty, numMismatches, numGapOpens, numGapExtend + 1);
            }

        }

        public static enum ExtensionType {
            MATCH, INSERTION, DELETION
        };

        public static enum Orientation {
            Forward("F"), Backward("B"), BidirectionalForward("BF"), BidirectionalBackward("BB");

            public final String symbol;

            private Orientation(String symbol) {
                this.symbol = symbol;
            }
        }

        public final ACGTSequence   read;

        public final Strand         strand;
        public final Orientation    orientation;
        public final ExtensionType  extensionType;
        public final int            cursor;
        public final Score          score;
        public final SuffixInterval si;

        protected Alignment(ACGTSequence read, Strand strand, Orientation orientation, ExtensionType extensionType,
                int cursor, Score score, SuffixInterval si) {
            this.read = read;
            this.strand = strand;
            this.orientation = orientation;
            this.extensionType = extensionType;
            this.cursor = cursor;
            this.score = score;
            this.si = si;
        }

        protected Alignment(Alignment other) {
            this.read = other.read;
            this.strand = other.strand;
            this.orientation = other.orientation;
            this.extensionType = other.extensionType;
            this.cursor = other.cursor;
            this.score = other.score;
            this.si = other.si;
        }

        public boolean gapIsAllowed(AlignmentScoreConfig config) {
            return score.numGapOpens < config.numGapOpenAllowed;
        }

        public boolean mismatchIsAllowed(AlignmentScoreConfig config) {
            return score.numMismatches < config.numMismatchesAllowed;
        }

        @Override
        public String toString() {
            return String.format("%s-%s(%d):%s/%s", orientation.symbol, extensionType.name().charAt(0), cursor,
                    score.toString(), strand.symbol);
        }

        public ACGT nextACGT() {
            return read.getACGT(cursor);
        }

        public ACGT getACGT(int cursor) {
            return read.getACGT(cursor);
        }

        public int getRemainingBases() {
            return Math.max(0, (int) read.textSize() - cursor);
        }

        public SuffixInterval suffixInterval() {
            return si;
        }

        public int getUpperBoundOfScore(AlignmentScoreConfig config) {
            return score.score + config.matchScore * getRemainingBases();
        }

        public boolean isFinished() {
            return cursor >= read.textSize();
        }

        public SuffixInterval nextSi(FMIndexOnGenome fmIndex, ACGT ch, SuffixInterval si) {
            switch (orientation) {
            case Forward:
            case BidirectionalForward:
                return fmIndex.forwardSearch(strand, ch, si);
            case Backward:
            case BidirectionalBackward:
                return fmIndex.backwardSearch(strand, ch, si);
            default:
                throw new IllegalStateException("cannot reach here");
            }
        }

        public abstract Alignment extend(ExtensionType type, int newCursor, Score newScore, SuffixInterval newSi);

        public Alignment extendWithMatch(AlignmentScoreConfig config, SuffixInterval newSi) {
            return extend(ExtensionType.MATCH, cursor + 1, score.extendWithMatch(config), newSi);
        }

        public Alignment extendWithMisMatch(AlignmentScoreConfig config, SuffixInterval newSi) {
            return extend(ExtensionType.MATCH, cursor + 1, score.extendWithMismatch(config), newSi);
        }

        public Alignment startInsertion(AlignmentScoreConfig config) {
            return extend(ExtensionType.INSERTION, cursor + 1, score.extendWithGapOpen(config), si);
        }

        public Alignment startDeletion(AlignmentScoreConfig config, SuffixInterval newSi) {
            return extend(ExtensionType.DELETION, cursor, score.extendWithGapOpen(config), newSi);
        }

        public Alignment extendInsertion(AlignmentScoreConfig config) {
            return extend(ExtensionType.INSERTION, cursor + 1, score.extendWithGapExtend(config), si);
        }

        public Alignment extendDeletion(AlignmentScoreConfig config, SuffixInterval newSi) {
            return extend(ExtensionType.DELETION, cursor, score.extendWithGapExtend(config), newSi);
        }

        public abstract Alignment newInstance(int newScore, int newCursor, SuffixInterval newSi);

        public abstract void accept(AlignmentVisitor visitor);

    }

    public static class ForwardAlignment extends Alignment
    {
        public static ForwardAlignment newInstance(ACGTSequence read, Strand strand, long N) {
            return new ForwardAlignment(read, strand, Orientation.Forward, ExtensionType.MATCH, 0, Score.initial(),
                    new SuffixInterval(0, N));
        }

        public ForwardAlignment(ACGTSequence read, Strand strand, Orientation orientation, ExtensionType extensionType,
                int cursor, Score score, SuffixInterval reverseSi) {
            super(read, strand, orientation, extensionType, cursor, score, reverseSi);
        }

        @Override
        public void accept(AlignmentVisitor visitor) {
            visitor.forwardAlignment(this);
        }

        @Override
        public SuffixInterval nextSi(FMIndexOnGenome fmIndex, ACGT ch, SuffixInterval si) {
            return fmIndex.forwardSearch(strand, ch, si);
        }

        @Override
        public Alignment extend(ExtensionType type, int newCursor, Score newScore, SuffixInterval newSi) {
            return new ForwardAlignment(read, strand, orientation, type, newCursor, newScore, newSi);
        }

        @Override
        public Alignment newInstance(int newScore, int newCursor, SuffixInterval newSi) {
            return new ForwardAlignment(read, strand, orientation, extensionType, newCursor, score.update(newScore),
                    newSi);
        }

    }

    public static class BackwardAlignment extends Alignment
    {
        public static BackwardAlignment newInstance(ACGTSequence read, Strand strand, long N) {
            return new BackwardAlignment(read, strand, Orientation.Backward, ExtensionType.MATCH, 0, Score.initial(),
                    new SuffixInterval(0, N));
        }

        public BackwardAlignment(ACGTSequence read, Strand strand, Orientation orientation,
                ExtensionType extensionType, int cursor, Score score, SuffixInterval forwardSi) {
            super(read, strand, orientation, extensionType, cursor, score, forwardSi);
        }

        @Override
        public ACGT getACGT(int cursor) {
            return read.getACGT((int) read.textSize() - 1 - cursor);
        }

        @Override
        public ACGT nextACGT() {
            return read.getACGT((int) read.textSize() - 1 - cursor);
        }

        @Override
        public void accept(AlignmentVisitor visitor) {
            visitor.backwardAlignment(this);
        }

        @Override
        public SuffixInterval nextSi(FMIndexOnGenome fmIndex, ACGT ch, SuffixInterval si) {
            return fmIndex.backwardSearch(strand, ch, si);
        }

        @Override
        public Alignment extend(ExtensionType type, int newCursor, Score newScore, SuffixInterval newSi) {
            return new BackwardAlignment(read, strand, orientation, type, newCursor, newScore, newSi);
        }

        @Override
        public Alignment newInstance(int newScore, int newCursor, SuffixInterval newSi) {
            return new BackwardAlignment(read, strand, orientation, extensionType, newCursor, score.update(newScore),
                    newSi);
        }

    }

    public static class BidirectionalForwardAlignment extends Alignment
    {
        public final SuffixInterval backwardSi;
        public final int            reverseCursor;

        public static BidirectionalForwardAlignment newInstance(ACGTSequence read, Strand strand, long N) {
            final int K = (int) read.textSize();
            int s2 = K / 3;
            return new BidirectionalForwardAlignment(read, strand, Orientation.BidirectionalForward,
                    ExtensionType.MATCH, s2, Score.initial(), new SuffixInterval(0, N), new SuffixInterval(0, N), s2);
        }

        protected BidirectionalForwardAlignment(ACGTSequence read, Strand strand, Orientation orientation,
                ExtensionType extensionType, int cursor, Score score, SuffixInterval forwardSi,
                SuffixInterval backwardSi, int reverseCursor) {
            super(read, strand, orientation, extensionType, cursor, score, forwardSi);
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
        public boolean isFinished() {
            return false;
        }

        @Override
        public void accept(AlignmentVisitor visitor) {
            visitor.bidirectionalForwardAlignment(this);
        }

        @Override
        public Alignment newInstance(int newScore, int newCursor, SuffixInterval newSi) {
            return new BidirectionalForwardAlignment(read, strand, orientation, extensionType, newCursor,
                    score.update(newScore), newSi, backwardSi, reverseCursor);
        }

        @Override
        public Alignment extend(ExtensionType type, int newCursor, Score newScore, SuffixInterval newSi) {
            // TODO fix reverse SI
            return new BidirectionalForwardAlignment(read, strand, orientation, type, newCursor, newScore, newSi,
                    backwardSi, newCursor);
        }

        @Override
        public SuffixInterval nextSi(FMIndexOnGenome fmIndex, ACGT ch, SuffixInterval si) {
            // TODO fix reverse Si
            SuffixInterval nextSi = fmIndex.forwardSearch(strand, ch, si);

            return nextSi;
        }
    }

}
