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
// ReadHit.java
// Since: 2011/08/25
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.record;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import org.utgenome.weaver.align.AlignmentScoreConfig;
import org.utgenome.weaver.align.CIGAR;
import org.utgenome.weaver.align.Strand;
import org.xerial.lens.SilkLens;

/**
 * Candidate position of read hit
 * 
 * @author leo
 * 
 */
public class ReadHit
{
    public final String chr;
    public final long   pos;
    public final int    matchLength;
    public final int    diff;
    public final Strand strand;
    public final CIGAR  cigar;
    public final int    numHits;
    public final int    qStart;
    public final int    qEnd;
    public ReadHit      nextSplit;

    public static ReadHit noHit(Strand strand) {
        return new ReadHit("", 0, 0, 0, 0, 0, strand, new CIGAR(), 0, null);
    }

    public ReadHit(String chr, long pos, int matchLength, int qStart, int qEnd, int diff, Strand strand, CIGAR cigar,
            int numHits, ReadHit nextSplit) {
        this.chr = chr;
        this.pos = pos;
        this.matchLength = matchLength;
        this.qStart = qStart;
        this.qEnd = qEnd;
        this.diff = diff;
        this.strand = strand;
        this.cigar = cigar;
        this.numHits = numHits;
        this.nextSplit = nextSplit;
    }

    /**
     * Match length including split alignment
     * 
     * @return
     */
    public int getTotalMatchLength() {
        if (nextSplit == null)
            return matchLength;
        else
            return matchLength + nextSplit.getTotalMatchLength();
    }

    public int getTotalDifferences() {
        if (nextSplit == null)
            return diff;
        else
            return diff + 1 + nextSplit.getTotalDifferences();
    }

    public boolean isUnique() {
        return numHits == 1;
    }

    public int getTotalScore(AlignmentScoreConfig config) {
        int fragmentLength = matchLength + diff;

        int score = matchLength * config.matchScore - diff * config.mismatchPenalty;
        if (nextSplit == null)
            return score;
        else
            return score - config.splitOpenPenalty + nextSplit.getTotalScore(config);
    }

    public String getAlignmentStateSingle() {
        return numHits > 0 ? ((numHits == 1) ? "U" : "R") : "N";
    }

    public String getAlignmentState() {
        String prefix = getAlignmentStateSingle();
        if (nextSplit == null)
            return prefix;
        else
            return prefix + nextSplit.getAlignmentState();
    }

    public String getAlignmentState(ReadHit head) {
        StringBuilder s = new StringBuilder();
        for (ReadHit c = head; c != null; c = c.nextSplit) {
            String state = c.getAlignmentStateSingle();
            if (c != this) {
                state = state.toLowerCase();
            }
            s.append(state);
        }
        return s.toString();
    }

    public String getCigarConcatenated() {
        if (nextSplit == null)
            return cigar.toString();
        else
            return cigar.toString() + "-" + nextSplit.getCigarConcatenated();
    }

    @Override
    public String toString() {
        return SilkLens.toSilk(this);
    }

    public ReadHit sortSplits() {
        if (this.nextSplit == null)
            return this;

        // Sort hits in the order of (chr, start) 
        ArrayList<ReadHit> r = new ArrayList<ReadHit>();
        for (ReadHit h = this; h != null; h = h.nextSplit) {
            r.add(h);
        }
        Collections.sort(r, new Comparator<ReadHit>() {
            @Override
            public int compare(ReadHit o1, ReadHit o2) {
                int diff = 0;
                if (o1.chr == null || o2.chr == null) {
                    diff = o1.qStart - o2.qStart;
                    if (!strand.isForward())
                        diff = -diff;
                }
                if (diff != 0)
                    return diff;

                diff = o1.chr.compareTo(o2.chr);
                if (diff != 0)
                    return diff;

                diff = (int) (o1.pos - o2.pos);
                return diff;
            }
        });

        // Fix split chain
        for (int i = 0; i < r.size() - 1; ++i) {
            r.get(i).nextSplit = r.get(i + 1);
        }
        r.get(r.size() - 1).nextSplit = null;
        return r.get(0);
    }

}
