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
// Genome Weaver Project
//
// SmithWatermanAlignment.java
// Since: Feb 22, 2010
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import org.utgenome.format.fasta.GenomeSequence;
import org.utgenome.weaver.align.CIGAR.Type;
import org.xerial.util.log.Logger;

/**
 * Smith-Waterman based aligner with Affine-gap score. This class implements
 * Smith-Waterman-Gotoh's local alignment algorithm. Standard/banded alignment
 * strategies are provided.
 * 
 * @author leo
 * 
 */
public class SmithWatermanAligner
{

    public static class Alignment
    {
        public final CIGAR          cigar;
        public final int            score;
        public final int            numMismatches;
        public final int            pos;          // 0-based leftmost position of the clipped sequence
        public final GenomeSequence rseq;         // reference sequence
        public final GenomeSequence qseq;         // query sequence

        public Alignment(CIGAR cigar, int score, int numMisamatches, GenomeSequence rseq, int pos, GenomeSequence qseq) {
            this.cigar = cigar;
            this.score = score;
            this.numMismatches = numMisamatches;
            this.pos = pos;
            this.rseq = rseq;
            this.qseq = qseq;
        }

        @Override
        public String toString() {
            return String.format("k:%d, cigar:%s, score:%d, pos:%d\n%s", numMismatches, cigar, score, pos,
                    detailedAlignment());
        }

        public String detailedAlignment() {
            if (rseq == null || qseq == null)
                return "";
            int max = Math.min(rseq.length(), qseq.length());

            int rCursor = 0;
            int rStart = pos;
            int qCursor = 0;

            StringBuilder r = new StringBuilder();
            StringBuilder s = new StringBuilder();
            StringBuilder q = new StringBuilder();

            if (cigar.size() > 0) {
                CIGAR.Element e = cigar.get(0);
                if (e.type == Type.SoftClip) {
                    int offset = pos - e.length;
                    rStart = offset;
                }
            }
            rCursor = Math.min(rCursor, rStart);
            for (; rCursor < rStart;) {
                r.append(rseq.charAt(rCursor++));
                s.append(" ");
                q.append(" ");
            }

            for (CIGAR.Element e : cigar.element()) {
                switch (e.type) {
                case Matches:
                case Mismatches:
                    for (int i = 0; i < e.length; ++i) {
                        ACGT rCh = ACGT.encode(getChar(rseq, rCursor++));
                        ACGT qCh = ACGT.encode(getChar(qseq, qCursor++));
                        r.append(rCh);
                        s.append(rCh.match(qCh) ? "|" : "X");
                        q.append(qCh);
                    }
                    break;
                case Insertions:
                    for (int i = 0; i < e.length; ++i) {
                        r.append("-");
                        s.append(" ");
                        q.append(getChar(qseq, qCursor++));
                    }
                    break;
                case Deletions:
                    for (int i = 0; i < e.length; ++i) {
                        r.append(getChar(rseq, rCursor++));
                        s.append(" ");
                        q.append("-");
                    }
                    break;
                case SoftClip:
                    for (int i = 0; i < e.length; ++i) {
                        r.append(getChar(rseq, rCursor++));
                        s.append(" ");
                        q.append(Character.toLowerCase(getChar(qseq, qCursor++)));
                    }
                    break;
                case Padding:
                case SkippedRegion:
                    for (int i = 0; i < e.length; ++i) {
                        r.append(getChar(rseq, rCursor++));
                        s.append(" ");
                        q.append(".");
                    }
                    break;
                case HardClip:
                    for (int i = 0; i < e.length; ++i) {
                        qCursor++;
                    }
                    break;
                }
            }

            return String.format("R: %s\n   %s\nQ: %s", r.toString(), s.toString(), q.toString());
        }

        private char getChar(GenomeSequence seq, int index) {
            return (index < seq.length() && index >= 0) ? seq.charAt(index) : ' ';
        }

    }

    public static enum Trace {
        NONE, DIAGONAL, DIAGONAL_MISMATCH, LEFT, UP
    };

    private static Logger              _logger = Logger.getLogger(SmithWatermanAligner.class);

    private final GenomeSequence       ref;
    private final GenomeSequence       query;
    private final AlignmentScoreConfig config;
    private final int                  N;                                                     // reference sequence length 
    private final int                  M;                                                     // query length
    private final int                  W;                                                     // band widdth
    private final int[][]              score;                                                 // DP matrix 
    // for gap-extension (Smith-Waterman Gotoh)
    private final int[][]              Li;
    private final int[][]              Ld;

    private int                        maxRow, maxCol, maxScore;                              // best score 

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

    protected void forwardDP(boolean bandedAlignment) {
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
    protected class DPScore
    {
        public final int M;
        public final int I;
        public final int D;
        private boolean  hasMatch;

        public DPScore(int row, int col) {
            ACGT r = ACGT.encode(ref.charAt(col - 1));
            ACGT q = ACGT.encode(query.charAt(row - 1));
            int scoreDiff;
            if (r.match(q)) {
                scoreDiff = config.matchScore;
                hasMatch = true;
            }
            else if (config.bssMode && r == ACGT.C && q == ACGT.T) {
                scoreDiff = -config.bssMismatchPenalty;
            }
            else {
                scoreDiff = -config.mismatchPenalty;
            }

            M = Math.max(score[row - 1][col - 1] + scoreDiff,
                    Math.max(Li[row - 1][col - 1] + scoreDiff, Ld[row - 1][col - 1] + scoreDiff));
            D = Math.max(score[row][col - 1] - config.gapOpenPenalty, Li[row][col - 1] - config.gapExtensionPenalty);
            I = Math.max(score[row - 1][col] - config.gapOpenPenalty, Ld[row - 1][col] - config.gapExtensionPenalty);
        }

        public boolean isMatch() {
            return hasMatch;
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
    protected Alignment traceback() {

        // trace back
        StringBuilder cigar = new StringBuilder();

        int leftMostPos = 0; // for seq1 

        // Append soft-clipped part in the query sequence
        for (int i = M - 1; i > maxRow; --i) {
            cigar.append("S");
        }

        int col = N - 1;
        int row = M - 1;

        // Append clipped sequences
        while (col > maxCol) {
            col--;
        }
        while (row > maxRow) {
            row--;
        }

        int diff = 0;

        // Trace back 
        traceback: for (col = maxCol, row = maxRow;;) {

            Trace path = Trace.NONE;
            boolean isMatch = false;
            // Recompute the score
            if (col >= 1 && row >= 1) {
                DPScore s = new DPScore(row, col);
                path = s.getPath();
                isMatch = s.isMatch();
            }

            switch (path) {
            case DIAGONAL:
                // match
                leftMostPos = col - 1;
                if (!isMatch) {
                    diff++;
                    cigar.append("X");
                }
                else {
                    cigar.append("M");
                }
                col--;
                row--;
                break;
            case LEFT:
                // insertion
                cigar.append("I");
                leftMostPos = col - 1;
                diff++;
                row--;
                break;
            case UP:
                // deletion
                cigar.append("D");
                diff++;
                col--;
                break;
            case NONE:
                while (col >= 1 || row >= 1) {
                    if (row >= 1) {
                        cigar.append("S");
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

        return new Alignment(new CIGAR(compactCigar.toString()), maxScore, diff, ref, leftMostPos, query);
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
    private static GenomeSequence wrap(String seq) {
        return new StringWrapper(seq);
    }

}
