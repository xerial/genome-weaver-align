/*--------------------------------------------------------------------------
 *  Copyright 2009 utgenome.org
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
// utgb-core Project
//
// SmithWatermanAlignment.java
// Since: Feb 22, 2010
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import org.utgenome.format.fasta.GenomeSequence;
import org.xerial.util.log.Logger;

/**
 * Smith-Waterman based aligner with Affine-gap score. Standard/banded alignment
 * is provided.
 * 
 * @author leo
 * 
 */
public class SmithWatermanAligner
{

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

        @Override
        public String toString() {
            return String.format("cigar:%s, score:%d, pos:%d\nrseq: %s\nqseq: %s", cigar, score, pos, rseq, qseq);
        }
    }

    private enum Trace {
        NONE, DIAGONAL, LEFT, UP
    };

    private static Logger              _logger = Logger.getLogger(SmithWatermanAligner.class);

    private final GenomeSequence       ref;
    private final GenomeSequence       query;
    private final AlignmentScoreConfig config;
    private final int                  N;
    private final int                  M;
    private final int                  W;
    private final int[][]              score;
    // for gap-extension (Smith-Waterman Gotoh)
    private final int[][]              Li;
    private final int[][]              Ld;

    private int                        maxRow, maxCol, maxScore;

    private SmithWatermanAligner(GenomeSequence ref, GenomeSequence query, AlignmentScoreConfig config) {
        this.ref = ref;
        this.query = query;
        this.config = config;
        this.N = ref.length() + 1;
        this.M = query.length() + 1;
        this.W = config.bandWidth;

        score = new int[M][N];
        Li = new int[M][N];
        Ld = new int[M][N];
    }

    public static Alignment standardAlign(GenomeSequence ref, GenomeSequence query) {
        return standardAlign(ref, query, new AlignmentScoreConfig());
    }

    public static Alignment standardAlign(GenomeSequence ref, GenomeSequence query, AlignmentScoreConfig config) {
        return align(ref, query, config, false);
    }

    public static Alignment bandedAlign(GenomeSequence ref, GenomeSequence query) {
        return bandedAlign(ref, query, new AlignmentScoreConfig());
    }

    public static Alignment bandedAlign(GenomeSequence ref, GenomeSequence query, AlignmentScoreConfig config) {
        return align(ref, query, config, true);
    }

    protected static Alignment align(GenomeSequence ref, GenomeSequence query, AlignmentScoreConfig config,
            boolean banded) {
        SmithWatermanAligner sw = new SmithWatermanAligner(ref, query, config);
        sw.forwardDP(banded);
        return sw.traceback();
    }

    private void forwardDP(boolean bandedAlignment) {
        // initialized the matrix
        final int MIN = Integer.MIN_VALUE / 2; // sufficiently small value 
        score[0][0] = 0;
        Li[0][0] = MIN;
        Ld[0][0] = MIN;
        for (int col = 1; col < N; ++col) {
            score[0][col] = 0; // Set 0 for local-alignment (Alignment can start from any position in the reference sequence)
            Li[0][col] = MIN;
            Ld[0][col] = -config.gapOpenPenalty - config.gapExtensionPenalty * (col - 1);
        }
        for (int row = 1; row < M; ++row) {
            score[row][0] = 0; // Setting this row to 0 allows clipped-alignment
            Li[row][0] = -config.gapOpenPenalty - config.gapExtensionPenalty * (row - 1);
            Ld[row][0] = MIN;
        }

        // dynamic programming
        int numProcessedCells = 0;
        int columnOffset = 1 - config.bandWidth / 2;
        for (int row = 1; row < M; ++row) {
            int colStart = bandedAlignment ? Math.max(1, columnOffset + row - 1) : 1;
            int columnMax = bandedAlignment ? (row == 1 ? N : Math.min(N, columnOffset + row - 1 + config.bandWidth))
                    : N;
            for (int col = colStart; col < columnMax; ++col) {
                DPScore s = new DPScore(row, col);
                if (!s.hasPositiveScore()) {
                    score[row][col] = 0;
                    Li[row][col] = 0;
                    Ld[row][col] = 0;
                    continue;
                }

                Li[row][col] = s.I;
                Ld[row][col] = s.D;
                score[row][col] = s.maxScore();

                // update max score
                if (score[row][col] > maxScore) {
                    maxRow = row;
                    maxCol = col;
                    maxScore = score[row][col];
                }
            }

        }

        //        if (_logger.isTraceEnabled())
        //            _logger.trace("N:%d, M:%d, NM:%d, num processed cells: %d", N, M, N * M, numProcessedCells);
    }

    /**
     * DP Score calculator
     * 
     * @author leo
     * 
     */
    private class DPScore
    {
        public final int M;
        public final int I;
        public final int D;

        public DPScore(int row, int col) {
            ACGT r = ACGT.encode(ref.charAt(col - 1));
            ACGT q = ACGT.encode(query.charAt(row - 1));
            int scoreDiff;
            if (r == q) {
                scoreDiff = config.matchScore;
            }
            else if (config.bssMode && r == ACGT.C && q == ACGT.T) {
                scoreDiff = -config.bssMismatchPenalty;
            }
            else {
                scoreDiff = -config.mismatchPenalty;
            }

            M = Math.max(score[row - 1][col - 1] + scoreDiff,
                    Math.max(Li[row - 1][col - 1] + scoreDiff, Ld[row - 1][col - 1] + scoreDiff));
            I = Math.max(score[row][col - 1] - config.gapOpenPenalty, Li[row][col - 1] - config.gapExtensionPenalty);
            D = Math.max(score[row - 1][col] - config.gapOpenPenalty, Ld[row - 1][col] - config.gapExtensionPenalty);
        }

        public int maxScore() {
            return Math.max(M, Math.max(I, D));
        }

        public boolean hasPositiveScore() {
            return M > 0 || I > 0 || D > 0;
        }

        public Trace getPath() {
            if (!hasPositiveScore())
                return Trace.NONE;

            if (M >= I) {
                if (M >= D)
                    return Trace.DIAGONAL;
                else
                    return Trace.UP;
            }
            else if (I >= D)
                return Trace.LEFT;
            else
                return Trace.UP;

        }
    }

    /**
     * Trace back the DP matrix, and obtain the alignment results
     * 
     * @return
     */
    private Alignment traceback() {

        // trace back
        StringBuilder cigar = new StringBuilder();
        StringBuilder a1 = new StringBuilder();
        StringBuilder a2 = new StringBuilder();

        int leftMostPos = 0; // for seq1 

        // Append soft-clipped part in the query sequence
        for (int i = M - 1; i > maxRow; --i) {
            cigar.append("S");
        }

        int col = N - 1;
        int row = M - 1;

        // Append clipped sequences
        while (col > maxCol) {
            a1.append(ref.charAt(col - 1));
            col--;
        }
        while (row > maxRow) {
            a2.append(Character.toLowerCase(query.charAt(row - 1)));
            row--;
        }

        // Trace back 
        traceback: for (col = maxCol, row = maxRow;;) {

            Trace path = Trace.NONE;
            // Recompute the score
            if (col >= 1 && row >= 1) {
                DPScore s = new DPScore(row, col);
                path = s.getPath();
            }

            switch (path) {
            case DIAGONAL:
                // match
                cigar.append("M");
                a1.append(ref.charAt(col - 1));
                a2.append(query.charAt(row - 1));
                leftMostPos = col - 1;
                col--;
                row--;
                break;
            case LEFT:
                // insertion
                cigar.append("I");
                a1.append("-");
                a2.append(query.charAt(row - 1));
                leftMostPos = col - 1;
                col--;
                break;
            case UP:
                cigar.append("D");
                a1.append(ref.charAt(col - 1));
                a2.append("-");
                row--;
                break;
            case NONE:
                while (col >= 1 || row >= 1) {
                    if (row >= 1) {
                        cigar.append("S");
                        a1.append(col >= 1 ? ref.charAt(col - 1) : ' ');
                        a2.append(Character.toLowerCase(query.charAt(row - 1)));
                    }
                    else {
                        a1.append(col >= 1 ? ref.charAt(col - 1) : ' ');
                        a2.append(' ');
                    }
                    col--;
                    row--;
                }

                break traceback; // exit the loop
            }
        }

        // create cigar string
        String cigarStr = cigar.reverse().toString();
        char prev = cigarStr.charAt(0);
        int count = 1;
        StringBuilder compactCigar = new StringBuilder();
        for (int i = 1; i < cigarStr.length(); ++i) {
            char c = cigarStr.charAt(i);
            if (prev == c) {
                count++;
            }
            else {
                compactCigar.append(Integer.toString(count));
                compactCigar.append(prev);

                prev = c;
                count = 1;
            }
        }
        if (count > 0) {
            compactCigar.append(Integer.toString(count));
            compactCigar.append(prev);
        }

        return new Alignment(compactCigar.toString(), maxScore, a1.reverse().toString(), leftMostPos, a2.reverse()
                .toString());
    }

    public static class StringWrapper implements GenomeSequence
    {

        public final String seq;

        public StringWrapper(String seq) {
            this.seq = seq;
        }

        public char charAt(int index) {
            return seq.charAt(index);
        }

        public int length() {
            return seq.length();
        }

    }

    /**
     * Wrap the input string as {@link GenomeSequence}
     * 
     * @param seq
     * @return
     */
    public static GenomeSequence wrap(String seq) {
        return new StringWrapper(seq);
    }

}
