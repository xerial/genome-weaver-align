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
// AlignmentState.java
// Since: 2011/02/17
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.util.ArrayList;
import java.util.List;

import org.utgenome.gwt.utgb.client.bio.CIGAR;
import org.utgenome.weaver.align.BWAlign.AlignmentScoreConfig;

public class Alignment implements Comparable<Alignment>
{
    public static enum IndelState {
        NORMAL, INSERTION, DELETION
    };

    public static class CommonInfo
    {
        public final IUPACSequence query;

        public CommonInfo(IUPACSequence query) {
            this.query = query;
        }
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

    public final CommonInfo           common;
    public final Strand               strand;
    public final int                  wordIndex;
    public final Alignment.IndelState indel;
    public final SuffixInterval       suffixInterval;
    public final int                  numMismatches;
    public final int                  alignmentScore;
    public final MismatchPosition     mismatchPosition;
    public final List<Gap>            gapPosition;

    private Alignment(CommonInfo common, Strand strand, int wordIndex, Alignment.IndelState indel,
            SuffixInterval suffixInterval, int numMismatches, int alignmentScore, MismatchPosition mismatchPosition,
            List<Gap> gapPosition) {
        this.common = common;
        this.strand = strand;
        this.wordIndex = wordIndex;
        this.indel = indel;
        this.suffixInterval = suffixInterval;
        this.numMismatches = numMismatches;
        this.alignmentScore = alignmentScore;
        this.mismatchPosition = mismatchPosition;
        this.gapPosition = gapPosition;
    }

    public CIGAR cigar() {
        int cursor = 0;
        CIGAR cigar = new CIGAR();
        if (gapPosition != null) {
            for (Gap gap : gapPosition) {
                int matchLen = gap.pos - cursor;
                if (matchLen > 0)
                    cigar.add(matchLen, CIGAR.Type.Matches);
                switch (gap.getType()) {
                case DELETION:
                    cigar.add(gap.len, CIGAR.Type.Deletions);
                    break;
                case INSERTION:
                    cigar.add(gap.len, CIGAR.Type.Insertions);
                    break;
                }
                cursor = gap.pos;
            }
        }

        int remainingLen = (int) (common.query.textSize() - cursor);
        if (remainingLen > 0)
            cigar.add(remainingLen, CIGAR.Type.Matches);

        return cigar;
    }

    public Alignment extendWithMatch(SuffixInterval next, AlignmentScoreConfig config) {
        return new Alignment(common, strand, this.wordIndex + 1, IndelState.NORMAL, next, numMismatches, alignmentScore
                + config.matchScore, mismatchPosition, gapPosition);
    }

    public Alignment extendWithMisMatch(SuffixInterval next, AlignmentScoreConfig config) {
        MismatchPosition newMM;
        if (mismatchPosition == null)
            newMM = MismatchPosition.oneBitInstance((int) (common.query.textSize() - 1), wordIndex + 1);
        else
            newMM = mismatchPosition.copyAndSet(wordIndex + 1);

        return new Alignment(common, strand, this.wordIndex + 1, IndelState.NORMAL, next, numMismatches + 1,
                alignmentScore - config.mismatchPenalty, newMM, gapPosition);
    }

    public Alignment extendWithDeletion(AlignmentScoreConfig config) {

        ArrayList<Gap> newGapPosition = new ArrayList<Gap>();

        int newScore = alignmentScore;
        if (indel == IndelState.DELETION) {
            assert (gapPosition != null);
            newScore -= config.gapExtentionPenalty;
            int i = 0;
            for (; i < gapPosition.size() - 1; ++i) {
                newGapPosition.add(gapPosition.get(i));
            }
            newGapPosition.add(gapPosition.get(i).extendOne());
        }
        else {
            newScore -= config.gapOpenPenalty;
            // update gap
            newGapPosition.add(new Deletion(wordIndex + 1, 1));
        }
        return new Alignment(common, strand, this.wordIndex, IndelState.DELETION, suffixInterval, numMismatches + 1,
                newScore, mismatchPosition, newGapPosition);
    }

    public Alignment extendWithInsertion(AlignmentScoreConfig config) {
        ArrayList<Gap> newGapPosition = new ArrayList<Gap>();

        int newScore = alignmentScore;
        if (indel == IndelState.INSERTION) {
            newScore -= config.gapExtentionPenalty;

            int i = 0;
            for (; i < gapPosition.size() - 1; ++i) {
                newGapPosition.add(gapPosition.get(i));
            }
            newGapPosition.add(gapPosition.get(i).extendOne());
        }
        else {
            newScore -= config.gapOpenPenalty;
            newGapPosition.add(new Insertion(wordIndex + 1, 1));
        }
        return new Alignment(common, strand, this.wordIndex + 1, IndelState.INSERTION, suffixInterval,
                numMismatches + 1, alignmentScore - config.gapExtentionPenalty, mismatchPosition, newGapPosition);
    }

    public static Alignment initialState(IUPACSequence seq, Strand strand, long textSize) {
        return new Alignment(new CommonInfo(seq), strand, 0, IndelState.NORMAL, new SuffixInterval(0, textSize - 1), 0,
                0, null, null);
    }

    @Override
    public int compareTo(Alignment o) {
        // Ascending order of the score
        return o.alignmentScore - this.alignmentScore;
    }

}
