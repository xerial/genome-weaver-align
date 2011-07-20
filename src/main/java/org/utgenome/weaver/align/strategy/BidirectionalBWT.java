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
import java.util.PriorityQueue;

import org.utgenome.weaver.align.ACGT;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.AlignmentSA;
import org.utgenome.weaver.align.AlignmentScoreConfig;
import org.utgenome.weaver.align.BitVector;
import org.utgenome.weaver.align.FMIndexOnGenome;
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

    private PriorityQueue<Alignment>   alignmentQueue = new PriorityQueue<Alignment>();

    public BidirectionalBWT(FMIndexOnGenome fmIndex, Reporter reporter) {
        this.fmIndex = fmIndex;
        this.reporter = reporter;
        this.N = fmIndex.textSize();
        this.config = new AlignmentScoreConfig();
    }

    public static class Alignment
    {
        public ACGTSequence   read;

        public SuffixInterval forward;
        public SuffixInterval backward;

        public int            cursorLeft;
        public int            cursorRight;
        public Strand         strand;

        public Alignment(ACGTSequence read, SuffixInterval forward, SuffixInterval backward, int cursorLeft,
                int cursorRight, Strand direction) {
            this.read = read;
            this.forward = forward;
            this.backward = backward;
            this.cursorLeft = cursorLeft;
            this.cursorRight = cursorRight;
            this.strand = direction;
        }

        public static Alignment startFromLeft(ACGTSequence read, Strand strand) {
            final int N = (int) read.textSize();
            return new Alignment(read, new SuffixInterval(0, N - 1), new SuffixInterval(0, N - 1), 0, 0, strand);
        }

        public static Alignment startFromMiddle(ACGTSequence read, Strand strand) {
            final int N = (int) read.textSize();
            return new Alignment(read, new SuffixInterval(0, N - 1), new SuffixInterval(0, N - 1), N / 2, N / 2, strand);
        }

        public static Alignment startFromRight(ACGTSequence read, Strand strand) {
            final int N = (int) read.textSize();
            return new Alignment(read, new SuffixInterval(0, N - 1), new SuffixInterval(0, N - 1), N - 1, N - 1, strand);
        }
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

        public QuickScanResult(SuffixInterval si, BitVector breakPoint, int numMismatches) {
            this.si = si;
            this.breakPoint = breakPoint;
            this.numMismatches = numMismatches;
        }
    }

    public QuickScanResult quickScan(ACGTSequence query, Strand direction) {
        int qLen = (int) query.textSize();
        int numMismatches = 0;
        BitVector breakPoint = new BitVector(qLen);
        SuffixInterval si = new SuffixInterval(0, N - 1);
        for (int i = 0; i < qLen; ++i) {
            ACGT ch = query.getACGT(i);
            si = fmIndex.backwardSearch(direction, ch, si);
            if (!si.isValidRange()) {
                si = new SuffixInterval(0, N - 1);
                breakPoint.set(i, true);
                numMismatches++;
            }
        }
        return new QuickScanResult(si, breakPoint, numMismatches);
    }

    void addQueue(Alignment alignment) {

    }

    void report(AlignmentSA result) throws Exception {
        fmIndex.toGenomeCoordinate(result, new ObjectHandlerBase<AlignmentRecord>() {
            @Override
            public void handle(AlignmentRecord r) throws Exception {
                reporter.emit(r);
            }
        });
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

    public void align(RawRead r) throws Exception {

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

        {
            // Break the read sequence into three parts 
            int s1 = (int) qF.textSize() / 3;
            int s2 = (int) qF.textSize() / 3 * 2;

            int[] m = new int[3];
            m[0] = scanF.breakPoint.countOneBits(0L, s1);
            m[1] = scanF.breakPoint.countOneBits(s1, s2);
            m[2] = scanF.breakPoint.countOneBits(s2, qF.textSize());

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
                    if (m[i - 1] < m[i])
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

            align(qF, Strand.FORWARD, searchStart);
        }

    }

    public void align(ACGTSequence query, Strand direction, SearchStart searchStart) {

        switch (searchStart) {
        case FRONT:
            break;
        case MIDDLE:
            break;
        case TAIL:
            break;
        }

    }

}
