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
import java.util.ArrayList;

import org.utgenome.UTGBException;
import org.utgenome.gwt.utgb.client.bio.IUPAC;
import org.utgenome.gwt.utgb.client.bio.SAMReadFlag;
import org.utgenome.weaver.align.SequenceBoundary.PosOnGenome;
import org.utgenome.weaver.align.record.AlignmentRecord;
import org.utgenome.weaver.align.strategy.BWAStrategy;
import org.xerial.lens.SilkLens;
import org.xerial.util.ObjectHandler;
import org.xerial.util.ObjectHandlerBase;
import org.xerial.util.StringUtil;
import org.xerial.util.log.Logger;
import org.xerial.util.opt.Argument;
import org.xerial.util.opt.Option;

/**
 * Burrows-Wheeler aligner for IUPAC sequences
 * 
 * @author leo
 * 
 */
public class BWAlign extends GenomeWeaverCommand
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

    //    @Option(longName = "sam", description = "output in SAM format")
    //    public boolean outputSAM           = false;

    @Option(symbol = "N", description = "Num mismatches allowed. default=0")
    public int     numMismachesAllowed = 0;

    @Override
    public void execute(String[] args) throws Exception {

        if (query == null) {
            throw new UTGBException("no query is given");
        }

        query(fastaFilePrefix, query, new ObjectHandlerBase<AlignmentRecord>() {
            @Override
            public void handle(AlignmentRecord r) throws Exception {

                ArrayList<Object> rec = new ArrayList<Object>();
                rec.add(r.readName);
                int flag = 0;
                if (r.strand == Strand.REVERSE)
                    flag |= SAMReadFlag.FLAG_STRAND_OF_QUERY;

                rec.add(flag);
                rec.add(r.chr);
                rec.add(r.start);
                rec.add(r.score);
                rec.add(r.getCigar().toCIGARString());
                rec.add("*");
                rec.add(0);
                rec.add(0);
                rec.add(r.querySeq);
                rec.add("*");
                rec.add("NM:i:" + r.numMismatches);
                System.out.println(StringUtil.join(rec, "\t"));
            }
        });
    }

    public static class FMIndexOnGenome
    {
        public final FMIndex            fmIndexF;
        public final FMIndex            fmIndexR;
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

        public void outputSAMHeader() {
            System.out.print(index.toSAMHeader());
        }

        public void toGenomeCoordinate(String querySeq, Alignment result, ObjectHandler<AlignmentRecord> reporter)
                throws Exception {
            //            if (_logger.isTraceEnabled())
            _logger.info(SilkLens.toSilk("alignment", result));

            final long querySize = result.common.query.textSize();

            for (long i = result.suffixInterval.lowerBound; i <= result.suffixInterval.upperBound; ++i) {
                long pos = -1;
                switch (result.strand) {
                case FORWARD:
                    pos = saR.get(i, fmIndexR);
                    pos = (fmIndexR.textSize() - 1 - pos) - querySize;
                    break;
                case REVERSE:
                    pos = saF.get(i, fmIndexF);
                    break;
                }
                pos += 1;
                if (pos != -1) {
                    PosOnGenome p = index.translate(pos, result.strand);
                    AlignmentRecord rec = new AlignmentRecord();
                    rec.chr = p.chr;
                    rec.start = p.pos;
                    rec.strand = result.strand;
                    rec.score = result.alignmentScore;
                    rec.numMismatches = result.numMismatches;
                    rec.querySeq = result.strand == Strand.FORWARD ? querySeq : result.common.query.reverse()
                            .toString();
                    rec.readName = querySeq;
                    rec.end = p.pos + result.wordIndex;
                    rec.setCIGAR(result.cigar().toCIGARString());
                    reporter.handle(rec);
                }
            }

        }
    }

    public static void query(String fastaFilePrefix, final String query,
            final ObjectHandler<AlignmentRecord> resultHandler) throws Exception {

        final FMIndexOnGenome fmIndex = new FMIndexOnGenome(fastaFilePrefix);

        BWAlign aligner = new BWAlign();
        aligner.fastaFilePrefix = fastaFilePrefix;
        aligner.query = query;

        _logger.info("query sequence: " + query);

        fmIndex.outputSAMHeader();

        BWAStrategy aln = new BWAStrategy(fmIndex, new ObjectHandlerBase<Alignment>() {
            @Override
            public void handle(Alignment input) throws Exception {
                fmIndex.toGenomeCoordinate(query, input, resultHandler);
            }
        });
        resultHandler.init();
        aln.align(query);
        resultHandler.finish();

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

}
