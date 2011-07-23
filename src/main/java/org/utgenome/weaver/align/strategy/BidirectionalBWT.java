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
import org.xerial.lens.JSONLens;
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

    public QuickScanResult quickScan(ACGTSequence query, Strand strand) {
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
                    return o2.getUpperBoundOfScore(AlignmentQueue.this.config)
                            - o1.getUpperBoundOfScore(AlignmentQueue.this.config);
                }
            });
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

            if (e.getUpperBoundOfScore(config) < scoreLowerBound) {
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
        QuickScanResult scanF = quickScan(qF, Strand.FORWARD);
        if (scanF.numMismatches == 0) {
            // Found an exact match
            report(AlignmentSA.exactMatch(config, r.name(), qF, scanF.si, Strand.FORWARD));
            return;
        }

        // Find potential mismatch positions for reverse direction
        ACGTSequence qC = qF.complement();
        QuickScanResult scanR = quickScan(qC, Strand.REVERSE);
        if (scanR.numMismatches == 0) {
            // Found an exact match
            report(AlignmentSA.exactMatch(config, r.name(), qC, scanR.si, Strand.REVERSE));
            return;
        }

        AlignmentQueue alignmentQueue = new AlignmentQueue(config);
        // Set the initial search states
        alignmentQueue.add(prepareInitialAlignmentState(qF, scanF, Strand.FORWARD));
        alignmentQueue.add(prepareInitialAlignmentState(qC, scanR, Strand.REVERSE));

        // Search iteration
        while (!alignmentQueue.isEmpty()) {
            Alignment current = alignmentQueue.poll();

            if (current.isFinished()) {
                report(current);
            }

            // extend the match
            switch (current.orientation) {
            case Forward: {
                int cursor = current.forwardCursor;
                SuffixInterval si = current.forwardSi;
                while (cursor < current.read.textSize() && si.isValidRange()) {
                    ACGT nextBase = current.read.getACGT(cursor);
                    SuffixInterval nextSi = fmIndex.forwardSearch(current.strand, nextBase, si);
                    if (nextSi.isValidRange()) {
                        si = nextSi;
                    }
                    else {
                        // Search for the other bases
                        for(ACGT ch : ACGT.values()) {
                            if(ch != nextBase)  {
                                SuffixInterval next = fmIndex.forwardSearch(current.strand, ch, si);
                                if(next.is)
                            }
                        }
                        
                    }
                    cursor++;
                }
                break;
            }
            case Backward:
                break;
            case BidirectionalBackward: {
                int cursor = current.backwardCursor;
                SuffixInterval si = current.reverseSi;
                while (cursor >= 0 && si.isValidRange()) {
                    ACGT nextBase = current.read.getACGT(cursor);
                    SuffixInterval nextSi = fmIndex.backwardSearch(current.strand, nextBase, si);
                    if (nextSi.isValidRange()) {
                        si = nextSi;
                    }
                    cursor--;
                }
            }
                break;
            }

        }

    }

    public static abstract class Alignment
    {
        public final ACGTSequence read;
        public final Strand       strand;
        public final int          cursor;
        public final int          score;
        public final int          numMismatches;
        public final int          numGapOpens;
        public final int          numGapExtend;

        protected Alignment(ACGTSequence read, Strand strand, int cursor, int score, int numMismatches,
                int numGapOpens, int numGapExtend) {
            this.read = read;
            this.strand = strand;
            this.cursor = cursor;
            this.score = score;
            this.numMismatches = numMismatches;
            this.numGapOpens = numGapOpens;
            this.numGapExtend = numGapExtend;
        }

        @Override
        public String toString() {
            return JSONLens.toJSON(this);
        }

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

    public static class ForwardAlignment extends Alignment
    {
        public final SuffixInterval reverseSi;

        public static ForwardAlignment newInstance(ACGTSequence read, Strand strand, long N) {
            return new ForwardAlignment(read, strand, 0, 0, 0, 0, 0, new SuffixInterval(0, N));
        }

        protected ForwardAlignment(ACGTSequence read, Strand strand, int cursor, int score, int numMismatches,
                int numGapOpens, int numGapExtend, SuffixInterval reverseSi) {
            super(read, strand, cursor, score, numMismatches, numGapOpens, numGapExtend);
            this.reverseSi = reverseSi;
        }

        @Override
        public int getUpperBoundOfScore(AlignmentScoreConfig config) {
            int remainingBases = Math.max(0, (int) read.textSize() - cursor);
            return score + config.matchScore * remainingBases;
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
            return new BackwardAlignment(read, strand, (int) read.textSize() - 1, 0, 0, 0, 0, new SuffixInterval(0, N));
        }

        public BackwardAlignment(ACGTSequence read, Strand strand, int cursor, int score, int numMismatches,
                int numGapOpens, int numGapExtend, SuffixInterval forwardSi) {
            super(read, strand, cursor, score, numMismatches, numGapOpens, numGapExtend);
            this.forwardSi = forwardSi;
        }

        @Override
        public int getUpperBoundOfScore(AlignmentScoreConfig config) {
            int remainingBases = Math.max(0, cursor + 1);
            return score + config.matchScore * remainingBases;
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
            return new BidirectionalForwardAlignment(read, strand, s2, 0, 0, 0, 0, new SuffixInterval(0, N),
                    new SuffixInterval(0, N), s2);
        }

        protected BidirectionalForwardAlignment(ACGTSequence read, Strand strand, int cursor, int score,
                int numMismatches, int numGapOpens, int numGapExtend, SuffixInterval forwardSi,
                SuffixInterval backwardSi, int reverseCursor) {
            super(read, strand, cursor, score, numMismatches, numGapOpens, numGapExtend);
            this.forwardSi = forwardSi;
            this.backwardSi = backwardSi;
            this.reverseCursor = reverseCursor;
        }

        @Override
        public int getUpperBoundOfScore(AlignmentScoreConfig config) {
            int remainingLeft = (int) read.textSize() / 3;
            int remainingRight = (int) read.textSize() - cursor;
            int remainingBases = remainingLeft + remainingRight;
            return score + config.matchScore * remainingBases;
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
