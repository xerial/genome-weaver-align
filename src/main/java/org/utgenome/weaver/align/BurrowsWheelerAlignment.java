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

import java.io.File;

import org.utgenome.gwt.utgb.client.bio.IUPAC;
import org.utgenome.weaver.align.IUPACSequence.IUPACBinaryInfo;
import org.xerial.util.opt.Argument;
import org.xerial.util.opt.Command;
import org.xerial.util.opt.Option;

/**
 * Burrows-Wheeler aligner for IUPAC sequences
 * 
 * @author leo
 * 
 */
public class BurrowsWheelerAlignment implements Command
{

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

    @Option(symbol = "L", description = "Save 1/L for occurrence count table (default = 256)")
    private int    L = 256;

    @Override
    public void execute(String[] args) throws Exception {

        IUPACBinaryInfo index = IUPACBinaryInfo.loadSilk(IUPACBinaryInfo.getFileName(fastaFilePrefix));
        final int N = index.totalSize;
        final int K = IUPAC.values().length;

        // loading BWT sequences
        File bwtForwardFile = new File(fastaFilePrefix + ".bwt");
        File bwtReverseFile = new File(fastaFilePrefix + ".rbwt");
        IUPACSequence bwtF = new IUPACSequence(bwtForwardFile, N);
        IUPACSequence bwtR = new IUPACSequence(bwtReverseFile, N);

        // Compute the occurrence tables
        OccurrenceCountTable occF = new OccurrenceCountTable(bwtF, L);
        OccurrenceCountTable occR = new OccurrenceCountTable(bwtR, L);

        // Count the character frequencies 
        CharacterCount C = new CharacterCount(bwtF);

    }

    public static class SuffixInterval
    {
        public final int lowerBound;
        public final int upperBound;

        public SuffixInterval(int lowerBound, int upperBound) {

            this.lowerBound = lowerBound;
            this.upperBound = upperBound;
        }

    }

    public static class AlignmentResult
    {
        public SuffixInterval interval;
    }

    public static interface Reporter<T>
    {
        public void emit(T result);
    }

    public static class FMIndexAlign
    {

        private final IUPACSequence             bwt;
        private final OccurrenceCountTable      occ;
        private final CharacterCount            C;
        private final Reporter<AlignmentResult> out;

        public FMIndexAlign(IUPACSequence bwt, OccurrenceCountTable occ, CharacterCount C, Reporter<AlignmentResult> out) {
            this.bwt = bwt;
            this.occ = occ;
            this.C = C;
            this.out = out;
        }

        public void align(String seq) {
            align(seq, seq.length() - 1, 3, new SuffixInterval(0, bwt.size() - 1));
        }

        public void align(String seq, int cursor, int numMismatchesAllowed, SuffixInterval si) {

            if (numMismatchesAllowed < 0)
                return;

            if (cursor < 0) {
                AlignmentResult result = new AlignmentResult();
                result.interval = si;
                out.emit(result);
                return;
            }

            // Search for deletion
            align(seq, cursor - 1, numMismatchesAllowed - 1, si);
            for (IUPAC base : IUPAC.values()) {
                int lowerBound = C.get(base) + occ.getOcc(base, si.lowerBound - 1);
                int upperBound = C.get(base) + occ.getOcc(base, si.upperBound) - 1;
                if (lowerBound < upperBound) {
                    SuffixInterval next = new SuffixInterval(lowerBound, upperBound);
                    // Search for insertion
                    align(seq, cursor, numMismatchesAllowed - 1, next);
                    IUPAC currentBase = IUPAC.encode(seq.charAt(cursor));

                    if ((base.bitFlag & currentBase.bitFlag) != 0) {
                        // match
                        align(seq, cursor - 1, numMismatchesAllowed, next);
                    }
                    else {
                        // mismatch
                        align(seq, cursor - 1, numMismatchesAllowed - 1, next);
                    }
                }
            }

        }

    }

}
