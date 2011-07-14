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
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.SuffixInterval;
import org.utgenome.weaver.align.record.ReadSequence;
import org.xerial.util.BitVector;

/**
 * Alignment algorithm using Bi-directional BWT
 * 
 * @author leo
 * 
 */
public class BidirectionalBWT
{
    private final FMIndexOnGenome    fmIndex;

    private PriorityQueue<Alignment> alignmentQueue = new PriorityQueue<Alignment>();

    public BidirectionalBWT(FMIndexOnGenome fmIndex) {
        this.fmIndex = fmIndex;
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

    void addQueue(Alignment alignment) {

    }

    public void align(ReadSequence read) {

        long N = fmIndex.textSize();
        int qLen = read.seq.length();
        ACGTSequence qF = new ACGTSequence(read.seq);
        ACGTSequence qC = qF.complement();

        // Find potential mismatch location
        BitVector mismatchPosition = new BitVector(qLen);
        SuffixInterval si = new SuffixInterval(0, N - 1);
        for (int i = 0; i < qLen; ++i) {
            ACGT ch = ACGT.encode(read.seq.charAt(i));
            si = fmIndex.forwardSearch(ch, si);
            if (!si.isValidRange()) {
                si = new SuffixInterval(0, N - 1);
                mismatchPosition.set(i, true);
            }
        }

    }

}
