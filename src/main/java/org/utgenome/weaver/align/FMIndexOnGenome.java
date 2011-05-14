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
// FMIndexOnGenome.java
// Since: 2011/04/28
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.io.IOException;
import java.io.PrintWriter;

import org.utgenome.UTGBException;
import org.utgenome.gwt.utgb.client.bio.IUPAC;
import org.utgenome.weaver.align.SequenceBoundary.PosOnGenome;
import org.utgenome.weaver.align.record.AlignmentRecord;
import org.xerial.lens.SilkLens;
import org.xerial.util.ObjectHandler;
import org.xerial.util.log.Logger;

public class FMIndexOnGenome
{
    private static Logger           _logger = Logger.getLogger(FMIndexOnGenome.class);

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
        BWAlign._logger.info("Loading sparse suffix arrays");
        saF = SparseSuffixArray.loadFrom(forwardDB.sparseSuffixArray());
        saR = SparseSuffixArray.loadFrom(reverseDB.sparseSuffixArray());

        // Load Wavelet arrays
        BWAlign._logger.info("Loading a Wavelet array of the forward BWT");
        wvF = WaveletArray.loadFrom(forwardDB.bwtWavelet());
        BWAlign._logger.info("Loading a Wavelet array of the reverse BWT");
        wvR = WaveletArray.loadFrom(reverseDB.bwtWavelet());

        // Prepare FM-indexes
        fmIndexF = new FMIndex(wvF);
        fmIndexR = new FMIndex(wvR);
    }

    public SuffixInterval backwardSearch(Strand strand, IUPAC nextBase, SuffixInterval si) {
        FMIndex fm = strand == Strand.FORWARD ? fmIndexR : fmIndexF;
        return fm.backwardSearch(nextBase, si);
    }

    public void outputSAMHeader(PrintWriter out) {
        out.print(index.toSAMHeader());
    }

    public long toForwardSequenceIndex(long saIndex, Strand strand) {
        long pos = -1;
        switch (strand) {
        case FORWARD:
            pos = fmIndexR.textSize() - 1 - saR.get(saIndex, fmIndexR);
            break;
        case REVERSE:
            pos = saF.get(saIndex, fmIndexF);
            break;
        }
        return pos;
    }

    public void toGenomeCoordinate(AlignmentSA result, ObjectHandler<AlignmentRecord> reporter) throws Exception {
        if (_logger.isTraceEnabled())
            _logger.info(SilkLens.toSilk("alignment", result));

        final long querySize = result.common.query.textSize();

        for (long i = result.suffixInterval.lowerBound; i <= result.suffixInterval.upperBound; ++i) {
            long pos = toForwardSequenceIndex(i, result.strand);
            if (result.strand == Strand.FORWARD)
                pos -= querySize;
            pos += 1;
            if (pos != -1) {
                PosOnGenome p = index.translate(pos, result.strand);
                AlignmentRecord rec = new AlignmentRecord();
                rec.chr = p.chr;
                rec.start = p.pos;
                rec.strand = result.strand;
                rec.score = result.alignmentScore;
                rec.numMismatches = result.numMismatches;
                // workaround for Picard tools, which cannot accept base character other than ACGT 
                rec.querySeq = result.strand == Strand.FORWARD ? result.common.query.toACGTString()
                        : result.common.query.reverse().toACGTString();
                rec.readName = result.common.queryName;
                rec.end = p.pos + result.wordIndex;
                rec.setCIGAR(result.cigar().toCIGARString());
                reporter.handle(rec);
            }
        }

    }
}
