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

import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.PriorityQueue;

import org.utgenome.UTGBException;
import org.utgenome.gwt.utgb.client.bio.IUPAC;
import org.utgenome.weaver.align.SequenceBoundary.PosOnGenome;
import org.xerial.lens.SilkLens;
import org.xerial.util.ObjectHandler;
import org.xerial.util.ObjectHandlerBase;
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

        if (query == null) {
            throw new UTGBException("no query is given");
        }
        final FMIndexOnGenome fmIndex = new FMIndexOnGenome(fastaFilePrefix);

        query(fastaFilePrefix, query, new ObjectHandlerBase<PosOnGenome>() {
            @Override
            public void handle(PosOnGenome input) throws Exception {
                _logger.info(SilkLens.toSilk(input));
            }
        });
    }

    public static class FMIndexOnGenome
    {
        private final FMIndex           fmIndexF;
        private final FMIndex           fmIndexR;
        private final SparseSuffixArray saF;
        private final SparseSuffixArray saR;
        private final WaveletArray      wvF;
        private final WaveletArray      wvR;
        private final SequenceBoundary  index;

        private final long              N;
        private final int               K;

        public FMIndexOnGenome(String fastaFilePrefix) throws UTGBException, IOException {
            BWTFiles forwardDB = new BWTFiles(fastaFilePrefix, Strand.FORWARD);
            BWTFiles reverseDB = new BWTFiles(fastaFilePrefix, Strand.REVERSE);

            // Load the boundary information of the concatenated chr sequences 
            index = SequenceBoundary.loadSilk(forwardDB.pacIndex());
            N = index.totalSize;
            K = IUPAC.values().length;

            // Load sparse suffix arrays
            _logger.info("Loading sparse suffix arrays");
            saF = SparseSuffixArray.loadFrom(forwardDB.sparseSuffixArray());
            saR = SparseSuffixArray.loadFrom(reverseDB.sparseSuffixArray());

            // Load Wavelet arrays
            _logger.info("Loading a Wavelet array of the forward BWT");
            wvF = WaveletArray.loadFrom(forwardDB.bwtWavelet());
            _logger.info("Loading a Wavelet array of the reverse BWT");
            wvR = WaveletArray.loadFrom(reverseDB.bwtWavelet());

            // Prepare FM-indexes
            fmIndexF = new FMIndex(wvF);
            fmIndexR = new FMIndex(wvR);
        }

        public void toGenomeCoordinate(Alignment result, ObjectHandler<PosOnGenome> reporter) throws Exception {
            _logger.info(SilkLens.toSilk("alignment", result));
            reporter.init();
            for (long i = result.suffixInterval.lowerBound; i <= result.suffixInterval.upperBound; ++i) {
                long pos = -1;
                switch (result.strand) {
                case FORWARD:
                    pos = saR.get(i, fmIndexR);
                    pos = (fmIndexF.textSize() - 1 - pos) - result.common.query.length();
                    break;
                case REVERSE:
                    pos = saF.get(i, fmIndexF);
                    break;
                }
                if (pos != -1)
                    reporter.handle(index.translate(pos));
            }
            reporter.finish();
        }
    }

    public static void query(String fastaFilePrefix, String query, final ObjectHandler<PosOnGenome> resultHandler)
            throws Exception {

        final FMIndexOnGenome fmIndex = new FMIndexOnGenome(fastaFilePrefix);

        BWAlign aligner = new BWAlign();
        aligner.fastaFilePrefix = fastaFilePrefix;
        aligner.query = query;

        _logger.info("query sequence: " + query);
        FMIndexAlign aln = new FMIndexAlign(fmIndex, new ObjectHandlerBase<Alignment>() {
            @Override
            public void handle(Alignment input) throws Exception {
                fmIndex.toGenomeCoordinate(input, resultHandler);
            }
        });
        aln.align(query);

    }

    static class FowardReverseFMIndex
    {
        final FMIndex F;
        final FMIndex R;

        public FowardReverseFMIndex(FMIndex forward, FMIndex reverse) {
            this.F = forward;
            this.R = reverse;
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

    /**
     * Alignment process using FM-Index
     * 
     * @author leo
     * 
     */
    public static class FMIndexAlign
    {
        private final FMIndexOnGenome          fmIndex;
        private final ObjectHandler<Alignment> out;

        private final PriorityQueue<Alignment> alignmentQueue       = new PriorityQueue<Alignment>();
        private final int                      numMismatchesAllowed = 1;

        private final AlignmentScoreConfig     config               = new AlignmentScoreConfig();

        private final ArrayList<IUPAC>         lettersInGenome      = new ArrayList<IUPAC>();

        private int                            bestScore            = -1;

        public FMIndexAlign(FMIndexOnGenome fmIndex, ObjectHandler<Alignment> out) {
            this.fmIndex = fmIndex;
            this.out = out;

            CharacterCount C = fmIndex.fmIndexF.getCharacterCount();
            for (IUPAC base : IUPAC.values()) {
                if (base == IUPAC.None)
                    continue;

                if (C.getCount(base) > 0) {
                    lettersInGenome.add(base);
                }
            }

        }

        public static String complement(String seq) {
            StringWriter rev = new StringWriter(seq.length());
            for (int i = 0; i < seq.length(); ++i) {
                char ch = Character.toUpperCase(seq.charAt(i));
                switch (ch) {
                case 'A':
                    rev.append('T');
                    break;
                case 'C':
                    rev.append('G');
                    break;
                case 'G':
                    rev.append('C');
                    break;
                case 'T':
                    rev.append('A');
                    break;
                default:
                case 'N':
                    rev.append('N');
                    break;
                // TODO IUPAC sequences
                }
            }
            return rev.toString();
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

            alignmentQueue.add(Alignment.initialState(seq, Strand.FORWARD, fmIndex.fmIndexR.textSize()));
            alignmentQueue.add(Alignment.initialState(complement(seq), Strand.REVERSE, fmIndex.fmIndexF.textSize()));

            while (!alignmentQueue.isEmpty()) {

                Alignment current = alignmentQueue.poll();
                if (current.numMismatches > numMismatchesAllowed) {
                    continue;
                }

                if (current.wordIndex >= seq.length()) {
                    if (current.alignmentScore >= bestScore) {
                        bestScore = current.alignmentScore;
                        out.handle(current);
                    }
                    continue;
                }

                // Search for deletion
                alignmentQueue.add(current.extendWithDeletion(config));
                //align(seq, cursor - 1, numMismatchesAllowed - 1, si);
                IUPAC currentBase = IUPAC.encode(current.common.query.charAt(current.wordIndex));
                for (IUPAC nextBase : lettersInGenome) {
                    FMIndex fm = current.strand == Strand.FORWARD ? fmIndex.fmIndexR : fmIndex.fmIndexF;
                    SuffixInterval next = fm.backwardSearch(nextBase, current.suffixInterval);
                    if (next.isValidRange()) {
                        // Search for insertion
                        if (current.wordIndex > 0 && current.wordIndex < seq.length() - 2) {
                            alignmentQueue.add(current.extendWithInsertion(config));
                        }
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
