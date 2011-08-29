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

/**
 * Simple Smith-Waterman aligner
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

    public static enum Trace {
        NONE, DIAGONAL, LEFT, UP
    };

    public static class Config
    {
        public int     MATCH_SCORE             = 1;
        public int     MISMATCH_PENALTY        = 3;
        public int     GAPOPEN_PENALTY         = 5;
        /**
         * for BSS alignment
         */
        public boolean BSS_ALIGNMENT           = true;
        public int     BSS_CT_MISMATCH_PENALTY = 0;
    }

    private final Config config;

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

    public Alignment align(String seq1, String seq2) {
        return align(wrap(seq1), wrap(seq2));
    }

    private final GenomeSequence seq1;
    private final GenomeSequence seq2;

    private final int            N;
    private final int            M;

    // prepare the score matrix
    private final int            score[][];
    private final Trace          trace[][];

    // maximum score and its location 
    private int                  maxScore = 0;
    private int                  maxX     = 0;
    private int                  maxY     = 0;

    private SmithWatermanAligner(GenomeSequence seq1, GenomeSequence seq2) {
        this(new Config(), seq1, seq2);
    }

    private SmithWatermanAligner(Config config, GenomeSequence seq1, GenomeSequence seq2) {
        this.config = config;
        this.seq1 = seq1;
        this.seq2 = seq2;
        this.N = seq1.length() + 1;
        this.M = seq2.length() + 1;

        score = new int[N][M];
        trace = new Trace[N][M];

        score[0][0] = 0;
        trace[0][0] = Trace.NONE;
    }

    public static void dpOnly(GenomeSequence seq1, GenomeSequence seq2) {
        SmithWatermanAligner sw = new SmithWatermanAligner(seq1, seq2);
        sw.align();
    }

    public static Alignment align(GenomeSequence seq1, GenomeSequence seq2) {
        SmithWatermanAligner sw = new SmithWatermanAligner(seq1, seq2);
        sw.align();
        return sw.traceBack();
    }

    private void align() {
        // set the first row and column all 0 for the local alignment 
        for (int x = 1; x < N; ++x) {
            score[x][0] = 0;
            trace[x][0] = Trace.NONE;
        }
        for (int y = 1; y < M; ++y) {
            score[0][y] = 0;
            trace[0][y] = Trace.NONE;
        }

        // SW loop
        for (int x = 1; x < N; ++x) {
            for (int y = 1; y < M; ++y) {
                char c1 = Character.toLowerCase(seq1.charAt(x - 1));
                char c2 = Character.toLowerCase(seq2.charAt(y - 1));

                // match(S), insertion(I) to , deletion(D) from the seq1 scores
                int S, I, D;
                if (c1 == c2)
                    S = score[x - 1][y - 1] + config.MATCH_SCORE;
                else {
                    if (config.BSS_ALIGNMENT && (c1 == 'c' && c2 == 't')) {
                        S = score[x - 1][y - 1] - config.BSS_CT_MISMATCH_PENALTY;
                    }
                    else
                        S = score[x - 1][y - 1] - config.MISMATCH_PENALTY;
                }

                I = score[x][y - 1] - config.GAPOPEN_PENALTY;
                D = score[x - 1][y] - config.GAPOPEN_PENALTY;

                if (S <= 0 && I <= 0 && D <= 0) {
                    score[x][y] = 0;
                    trace[x][y] = Trace.NONE;
                    continue;
                }

                // choose the best score
                if (S >= I) {
                    if (S >= D) {
                        score[x][y] = S;
                        trace[x][y] = Trace.DIAGONAL;
                    }
                    else {
                        score[x][y] = D;
                        trace[x][y] = Trace.LEFT;
                    }
                }
                else if (I >= D) {
                    score[x][y] = I;
                    trace[x][y] = Trace.UP;
                }
                else {
                    score[x][y] = D;
                    trace[x][y] = Trace.LEFT;
                }

                // update max score
                if (score[x][y] > maxScore) {
                    maxX = x;
                    maxY = y;
                    maxScore = score[x][y];
                }
            }
        }
    }

    private Alignment traceBack() {

        // trace back
        StringBuilder cigar = new StringBuilder();
        StringBuilder a1 = new StringBuilder();
        StringBuilder a2 = new StringBuilder();

        int leftMostPos = 0; // for seq1 

        for (int i = M - 1; i > maxY; --i) {
            cigar.append("S");
        }
        boolean toContinue = true;

        int x = N - 1;
        int y = M - 1;

        while (x > maxX) {
            a1.append(seq1.charAt(x - 1));
            x--;
        }
        while (y > maxY) {
            a2.append(Character.toLowerCase(seq2.charAt(y - 1)));
            y--;
        }

        for (x = maxX, y = maxY; toContinue;) {
            switch (trace[x][y]) {
            case DIAGONAL:
                cigar.append("M");
                a1.append(seq1.charAt(x - 1));
                a2.append(seq2.charAt(y - 1));
                leftMostPos = x - 1;
                x--;
                y--;
                break;
            case LEFT:
                cigar.append("D");
                a1.append("-");
                a2.append(seq2.charAt(y - 1));
                leftMostPos = x - 1;
                x--;
                break;
            case UP:
                cigar.append("I");
                a1.append(seq1.charAt(x - 1));
                a2.append("-");
                y--;
                break;
            case NONE:
                toContinue = false;
                while (x >= 1 || y >= 1) {
                    if (y >= 1) {
                        cigar.append("S");
                        a1.append(x >= 1 ? seq1.charAt(x - 1) : ' ');
                        a2.append(Character.toLowerCase(seq2.charAt(y - 1)));
                    }
                    else {
                        a1.append(x >= 1 ? seq1.charAt(x - 1) : ' ');
                        a2.append(' ');
                    }
                    x--;
                    y--;
                }

                break;
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

}
