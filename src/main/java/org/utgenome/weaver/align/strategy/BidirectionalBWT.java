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

import java.util.PriorityQueue;

import org.utgenome.weaver.align.ACGT;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.AlignmentSA;
import org.utgenome.weaver.align.AlignmentScoreConfig;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.SuffixInterval;
import org.utgenome.weaver.align.record.AlignmentRecord;
import org.utgenome.weaver.align.record.RawRead;
import org.utgenome.weaver.align.record.ReadSequence;
import org.utgenome.weaver.parallel.Reporter;
import org.xerial.util.BitVector;
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
        public final BitVector      mismatchPosition;
        public final int            numMismatches;

        public QuickScanResult(SuffixInterval si, BitVector mismatchPosition, int numMismatches) {
            this.si = si;
            this.mismatchPosition = mismatchPosition;
            this.numMismatches = numMismatches;
        }
    }

    public QuickScanResult quickScan(ACGTSequence query, Strand direction) {
        int qLen = (int) query.textSize();
        int numMismatches = 0;
        BitVector mismatchPosition = new BitVector(qLen);
        SuffixInterval si = new SuffixInterval(0, N - 1);
        for (int i = 0; i < qLen; ++i) {
            ACGT ch = query.getACGT(i);
            si = fmIndex.backwardSearch(direction, ch, si);
            if (!si.isValidRange()) {
                si = new SuffixInterval(0, N - 1);
                mismatchPosition.set(i, true);
                numMismatches++;
            }
        }
        return new QuickScanResult(si, mismatchPosition, numMismatches);
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

    public void align(RawRead r) throws Exception {

        ReadSequence read = (ReadSequence) r;

        ACGTSequence qF = new ACGTSequence(read.seq);

        // Find potential mismatch positions for forward direction
        QuickScanResult scanF = quickScan(qF, Strand.FORWARD);
        if (scanF.numMismatches == 0) {
            // Found an exact match
            AlignmentSA result = AlignmentSA.exactMatch(config, r.name(), qF, scanF.si, Strand.FORWARD);
            report(result);
            return;
        }

        // Find potential mismatch positions for reverse direction
        ACGTSequence qC = qF.complement();
        QuickScanResult scanR = quickScan(qC, Strand.REVERSE);
        if (scanR.numMismatches == 0) {
            // Found an exact match
            AlignmentSA result = AlignmentSA.exactMatch(config, r.name(), qC, scanR.si, Strand.REVERSE);
            report(result);
            return;
        }

        if (scanF.numMismatches < scanR.numMismatches) {
            // Search 

        }
        else {

        }

    }

}
