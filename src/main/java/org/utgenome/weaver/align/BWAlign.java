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
// BurrowsWheelerAlignment.java
// Since: 2011/02/10
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.PriorityQueue;

import org.utgenome.gwt.utgb.client.bio.IUPAC;
import org.utgenome.weaver.align.SequenceBoundary.PosOnGenome;
import org.xerial.lens.SilkLens;
import org.xerial.util.log.Logger;
import org.xerial.util.opt.Argument;
import org.xerial.util.opt.Command;
import org.xerial.util.opt.Option;

/**
 * Burrows-Wheeler aligner for IUPAC sequences
 * 
 * @author leo
 * 
 */
public class BWAlign implements Command
{
    private static Logger _logger = Logger.getLogger(BWAlign.class);

    @Override
    public String name() {
        return "align";
    }

    @Override
    public String getOneLineDescription() {
        return "performs alignment";
    }

    @Override
    public Object getOptionHolder() {
        return this;
    }

    @Argument(index = 0)
    private String fastaFilePrefix;

    @Option(symbol = "L", description = "1/L of the Occ table will be precomputed. (default L = 256)")
    private int    L = 256;

    @Option(symbol = "q", description = "query sequence")
    private String query;

    @Override
    public void execute(String[] args) throws Exception {

        SequenceBoundary index = SequenceBoundary.loadSilk(SequenceBoundary.getFileName(fastaFilePrefix));
        final int N = index.totalSize;
        final int K = IUPAC.values().length;

        // Load BWT sequences
        File bwtForwardFile = new File(fastaFilePrefix + ".bwt");
        File bwtReverseFile = new File(fastaFilePrefix + ".rbwt");
        IUPACSequence bwtF = new IUPACSequence(bwtForwardFile, N);
        IUPACSequence bwtR = new IUPACSequence(bwtReverseFile, N);

        // Load the boundary information of the concatenated chr sequences 
        final SequenceBoundary boundary = SequenceBoundary.loadSilk(SequenceBoundary.getFileName(fastaFilePrefix));

        // Load sparse suffix arrays
        File sparseForwardSAFile = new File(fastaFilePrefix + ".sa");
        File sparseReverseSAFile = new File(fastaFilePrefix + ".rsa");
        final SparseSuffixArray saF = SparseSuffixArray.loadFrom(new BufferedInputStream(new FileInputStream(
                sparseForwardSAFile)));
        final SparseSuffixArray saR = SparseSuffixArray.loadFrom(new BufferedInputStream(new FileInputStream(
                sparseReverseSAFile)));

        // Compute the occurrence tables
        OccurrenceCountTable occF = new OccurrenceCountTable(bwtF, L);
        OccurrenceCountTable occR = new OccurrenceCountTable(bwtR, L);

        // Count the character frequencies 
        CharacterCount C = new CharacterCount(bwtF);

        if (query != null) {
            _logger.info("query sequence: " + query);
            final FMIndex fmIndex = new FMIndex(bwtF, occF, C);
            FMIndexAlign aln = new FMIndexAlign(fmIndex, new Reporter<AlignmentResult>() {
                @Override
                public void emit(AlignmentResult result) throws Exception {
                    _logger.info(SilkLens.toSilk("alignment", result));
                    for (int i = result.suffixInterval.lowerBound; i <= result.suffixInterval.upperBound; ++i) {
                        int pos = saF.get(i, fmIndex);
                        PosOnGenome loc = boundary.translate(pos);
                        _logger.info(SilkLens.toSilk("loc", loc));
                    }
                }
            });

            aln.align(query);
        }

    }

    public static String reverse(String query) {
        final int N = query.length();
        StringBuilder buf = new StringBuilder(N);
        for (int i = N - 1; i >= 0; --i) {
            buf.append(query.charAt(i));
        }
        return buf.toString();
    }

    public static class AlignmentResult
    {
        public SuffixInterval suffixInterval;
        public int            numMismatches = 0;
    }

    public static interface Reporter<T>
    {
        public void emit(T result) throws Exception;
    }

    public static class Gap
    {
        public int pos;
    }

    public static class AlignmentScoreConfig
    {
        public final int matchScore          = 1;
        public final int mismatchPenalty     = 3;
        public final int gapOpenPenalty      = 11;
        public final int gapExtentionPenalty = 4;
    }

    public static class AlignmentState implements Comparable<AlignmentState>
    {
        public int            wordIndex;
        public SuffixInterval suffixInterval;
        public int            numMismatches;
        public int            alignmentScore;
        public List<Integer>  mismatchPosition;
        public List<Integer>  gapPosition;

        private AlignmentState(int wordIndex, SuffixInterval suffixInterval, int numMismatches, int alignmentScore,
                List<Integer> mismatchPosition, List<Integer> gapPosition) {
            this.wordIndex = wordIndex;
            this.suffixInterval = suffixInterval;
            this.numMismatches = numMismatches;
            this.alignmentScore = alignmentScore;
            this.mismatchPosition = mismatchPosition;
            this.gapPosition = gapPosition;
        }

        public AlignmentState extendWithMatch(SuffixInterval next, AlignmentScoreConfig config) {
            return new AlignmentState(this.wordIndex - 1, next, numMismatches, alignmentScore + config.matchScore,
                    mismatchPosition, gapPosition);
        }

        public AlignmentState extendWithMisMatch(SuffixInterval next, AlignmentScoreConfig config) {
            ArrayList<Integer> newMismatchPosition = new ArrayList<Integer>();
            if (mismatchPosition != null)
                newMismatchPosition.addAll(mismatchPosition);
            newMismatchPosition.add(wordIndex + 1);
            return new AlignmentState(this.wordIndex - 1, next, numMismatches + 1, alignmentScore
                    - config.mismatchPenalty, newMismatchPosition, gapPosition);
        }

        public AlignmentState extendWithDeletion(AlignmentScoreConfig config) {
            ArrayList<Integer> newGapPosition = new ArrayList<Integer>();
            if (gapPosition != null)
                newGapPosition.addAll(gapPosition);
            newGapPosition.add(wordIndex + 1);
            return new AlignmentState(this.wordIndex - 1, suffixInterval, numMismatches + 1, alignmentScore
                    - config.gapExtentionPenalty, mismatchPosition, gapPosition);
        }

        public AlignmentState extendWithInsertion(AlignmentScoreConfig config) {
            ArrayList<Integer> newGapPosition = new ArrayList<Integer>();
            if (gapPosition != null)
                newGapPosition.addAll(gapPosition);
            // TODO distinguish indels
            newGapPosition.add(wordIndex + 1);
            return new AlignmentState(this.wordIndex - 1, suffixInterval, numMismatches + 1, alignmentScore
                    - config.gapExtentionPenalty, mismatchPosition, gapPosition);
        }

        public static AlignmentState initialState(String seq, FMIndex fmIndex) {
            return new AlignmentState(seq.length() - 1, new SuffixInterval(0, fmIndex.textSize() - 1), 0, 0, null, null);
        }

        @Override
        public int compareTo(AlignmentState o) {
            // Ascending order of the score
            return o.alignmentScore - this.alignmentScore;
        }

    }

    public static class FMIndexAlign
    {
        private final FMIndex                       fmIndex;
        private final Reporter<AlignmentResult>     out;

        private final PriorityQueue<AlignmentState> alignmentQueue       = new PriorityQueue<AlignmentState>();
        private final int                           numMismatchesAllowed = 0;

        private final AlignmentScoreConfig          config               = new AlignmentScoreConfig();

        public FMIndexAlign(FMIndex fmIndex, Reporter<AlignmentResult> out) {
            this.fmIndex = fmIndex;
            this.out = out;
        }

        public void align(String seq) throws Exception {

            alignmentQueue.add(AlignmentState.initialState(seq, fmIndex));
            alignSeq(seq);
        }

        /**
         * 
         * @param seq
         * @param cursor
         * @param numMismatchesAllowed
         * @param si
         * @throws Exception
         */
        public void alignSeq(String seq) throws Exception {

            while (!alignmentQueue.isEmpty()) {

                AlignmentState current = alignmentQueue.poll();
                if (current.numMismatches > numMismatchesAllowed) {
                    continue;
                }

                if (current.wordIndex < 0) {
                    AlignmentResult result = new AlignmentResult();
                    result.suffixInterval = current.suffixInterval;
                    result.numMismatches = current.numMismatches;
                    out.emit(result);
                    continue;
                }

                // Search for deletion
                alignmentQueue.add(current.extendWithDeletion(config));
                //align(seq, cursor - 1, numMismatchesAllowed - 1, si);
                IUPAC currentBase = IUPAC.encode(seq.charAt(current.wordIndex));
                for (IUPAC nextBase : new IUPAC[] { IUPAC.A, IUPAC.C, IUPAC.G, IUPAC.T }) {
                    SuffixInterval next = fmIndex.backwardSearch(nextBase, current.suffixInterval);
                    if (next.isValidRange()) {
                        // Search for insertion
                        alignmentQueue.add(current.extendWithInsertion(config));
                        if ((nextBase.bitFlag & currentBase.bitFlag) != 0) {
                            // match
                            alignmentQueue.add(current.extendWithMatch(next, config));
                        }
                        else {
                            // mismatch
                            alignmentQueue.add(current.extendWithMisMatch(next, config));
                        }
                    }
                }
            }

        }

    }

}
