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

    @Option(symbol = "q", description = "query sequence")
    private String query;

    @Override
    public void execute(String[] args) throws Exception {

        BWTFiles forwardDB = new BWTFiles(fastaFilePrefix, Strand.FORWARD);
        BWTFiles reverseDB = new BWTFiles(fastaFilePrefix, Strand.REVERSE);

        // Load the boundary information of the concatenated chr sequences 
        final SequenceBoundary index = SequenceBoundary.loadSilk(forwardDB.pacIndex());
        final long N = index.totalSize;
        final int K = IUPAC.values().length;

        // Load sparse suffix arrays
        _logger.info("Loading sparse suffix arrays");
        final SparseSuffixArray saF = SparseSuffixArray.loadFrom(forwardDB.sparseSuffixArray());
        final SparseSuffixArray saR = SparseSuffixArray.loadFrom(reverseDB.sparseSuffixArray());

        // Load Wavelet arrays
        _logger.info("Loading a Wavelet array of the forward BWT");
        WaveletArray wvF = WaveletArray.loadFrom(forwardDB.bwtWavelet());
        _logger.info("Loading a Wavelet array of the reverse BWT");
        WaveletArray wvR = WaveletArray.loadFrom(reverseDB.bwtWavelet());

        final FMIndex fmIndex = new FMIndex(wvR);

        if (query != null) {
            _logger.info("query sequence: " + query);

            FMIndexAlign aln = new FMIndexAlign(fmIndex, new Reporter<AlignmentState>() {
                @Override
                public void emit(AlignmentState result) throws Exception {
                    _logger.info(SilkLens.toSilk("alignment", result));
                    for (long i = result.suffixInterval.lowerBound; i <= result.suffixInterval.upperBound; ++i) {
                        long pos = saR.get(i, fmIndex);
                        if (result.strand == Strand.FORWARD) {
                            pos = (fmIndex.textSize() - 1 - pos) - query.length();
                        }
                        PosOnGenome loc = index.translate(pos);
                        System.out.println(SilkLens.toSilk("loc", loc));
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

    public static interface Reporter<T>
    {
        public void emit(T result) throws Exception;
    }

    public static class AlignmentScoreConfig
    {
        public final int matchScore          = 1;
        public final int mismatchPenalty     = 3;
        public final int gapOpenPenalty      = 11;
        public final int gapExtentionPenalty = 4;
    }

    public static abstract class Gap
    {
        public static enum Type {
            INSERTION, DELETION
        }

        public final int pos;
        public final int len;

        public Gap(int pos, int len) {
            this.pos = pos;
            this.len = len;
        }

        public abstract Gap extendOne();

        public abstract Type getType();
    }

    public static class Insertion extends Gap
    {

        public Insertion(int pos, int len) {
            super(pos, len);
        }

        @Override
        public Gap extendOne() {
            return new Insertion(pos, len + 1);
        }

        @Override
        public Type getType() {
            return Type.INSERTION;
        }
    }

    public static class Deletion extends Gap
    {

        public Deletion(int pos, int len) {
            super(pos, len);
        }

        @Override
        public Gap extendOne() {
            return new Deletion(pos, len + 1);
        }

        @Override
        public Type getType() {
            return Type.INSERTION;
        }
    }

    public static class FMIndexAlign
    {
        private final FMIndex                       fmIndex;
        private final Reporter<AlignmentState>      out;

        private final PriorityQueue<AlignmentState> alignmentQueue       = new PriorityQueue<AlignmentState>();
        private final int                           numMismatchesAllowed = 0;

        private final AlignmentScoreConfig          config               = new AlignmentScoreConfig();

        public FMIndexAlign(FMIndex fmIndex, Reporter<AlignmentState> out) {
            this.fmIndex = fmIndex;
            this.out = out;
        }

        /**
         * 
         * @param seq
         * @param cursor
         * @param numMismatchesAllowed
         * @param si
         * @throws Exception
         */
        public void align(String seq) throws Exception {

            alignmentQueue.add(AlignmentState.initialState(seq, Strand.FORWARD, fmIndex));

            while (!alignmentQueue.isEmpty()) {

                AlignmentState current = alignmentQueue.poll();
                if (current.numMismatches > numMismatchesAllowed) {
                    continue;
                }

                if (current.wordIndex >= seq.length()) {
                    out.emit(current);
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
