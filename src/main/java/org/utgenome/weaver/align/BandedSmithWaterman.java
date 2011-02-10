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
// BandedSmithWaterman.java
// Since: 2011/01/31
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;


/**
 * Banded Smith-Waterman Alignment
 * 
 * @author leo
 * 
 */
public class BandedSmithWaterman
{
    private final int     ROW_SIZE;
    private final int     COL_SIZE;
    private final int     BAND_WIDTH;
    private final int[][] bandedScoringMatrix;

    public BandedSmithWaterman(int rowSize, int colSize, int bandWidth) {
        this.ROW_SIZE = rowSize;
        this.COL_SIZE = colSize;
        this.BAND_WIDTH = bandWidth;

        bandedScoringMatrix = new int[COL_SIZE][bandWidth];
    }

    public static class Alignment
    {
        public final String cigar;
        public final int    score;
        public final int    pos;  // 1-based leftmost position of the clipped sequence
        public final String rseq; // reference sequence
        public final String qseq; // query sequence

        public Alignment(String cigar, int score, String rseq, int pos, String qseq) {
            this.cigar = cigar;
            this.score = score;
            this.pos = pos;
            this.rseq = rseq;
            this.qseq = qseq;
        }
    }

    //    public Alignment align(GenomeSequence seq1, GenomeSequence seq2) {
    //
    //    }

}
