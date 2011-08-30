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
// AlignmentRecord.java
// Since: 2011/04/25
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.record;

import java.util.ArrayList;

import org.utgenome.gwt.utgb.client.UTGBClientException;
import org.utgenome.gwt.utgb.client.bio.CIGAR;
import org.utgenome.gwt.utgb.client.bio.CIGAR.Type;
import org.utgenome.gwt.utgb.client.bio.SAMReadFlag;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.Strand;
import org.xerial.util.StringUtil;

/**
 * Holds an alignment result in a form that can be converted into SAM format
 * 
 * @author leo
 * 
 */
public class AlignmentRecord
{
    public String          readName;
    public String          chr;
    public Strand          strand;
    public int             start;
    public int             end;
    public int             numMismatches = 0;
    private CIGAR          cigar;
    public String          querySeq;
    public String          qual;
    public int             score;
    public int             numBestHits;
    public AlignmentRecord split         = null;

    public AlignmentRecord() {

    }

    public AlignmentRecord(String readName, String chr, Strand strand, int start, int end, int numMismatches,
            CIGAR cigar, String querySeq, String qual, int score, int numBestHists, AlignmentRecord split) {
        this.readName = readName;
        this.chr = chr;
        this.strand = strand;
        this.start = start;
        this.end = end;
        this.numMismatches = numMismatches;
        this.cigar = cigar;
        this.querySeq = querySeq;
        this.qual = qual;
        this.score = score;
        this.numBestHits = numBestHists;
        this.split = split;
    }

    public CIGAR getCigar() {
        if (cigar != null)
            return cigar;
        else
            return new CIGAR();
    }

    public void setCIGAR(String cigarStr) throws UTGBClientException {
        this.cigar = new CIGAR(cigarStr);
    }

    public String toSAMLine() {
        ArrayList<Object> column = new ArrayList<Object>();
        column.add(readName);
        int flag = 0;
        if (strand == Strand.REVERSE)
            flag |= SAMReadFlag.FLAG_STRAND_OF_QUERY;

        column.add(flag);
        column.add(chr);
        column.add(start);
        column.add(score);
        column.add(getCigar());
        column.add("*"); // pair chr
        column.add(0); // pair start
        column.add(0); // insert size
        column.add(querySeq);
        column.add(qual == null ? "*" : qual); // quality value
        if (numMismatches >= 0)
            column.add("NM:i:" + numMismatches);
        String line = StringUtil.join(column, "\t");
        if (split != null)
            line += "\n" + split.toSAMLine();

        return line;
    }

    private static String reverse(String s) {
        if (s == null)
            return null;
        StringBuilder out = new StringBuilder(s.length());
        for (int i = s.length() - 1; i >= 0; --i)
            out.append(s.charAt(i));
        return out.toString();
    }

    public static AlignmentRecord convert(ReadHit hit, Read read, int numOtherBestHits) throws UTGBClientException {

        final int numHits = hit.numHits + numOtherBestHits;
        ACGTSequence query = read.getRead(0);
        final int m = (int) query.textSize();

        String qual = read.getQual(0);
        if (qual == null)
            qual = "";
        if (!hit.strand.isForward()) {
            query = query.reverseComplement();
            qual = reverse(qual);
        }

        if (hit.nextSplit == null) {
            // TODO alignment score
            CIGAR cigar = new CIGAR(hit.cigar);
            //cigar.add(hit.matchLength, Type.Matches);
            AlignmentRecord rec = new AlignmentRecord(read.name(), hit.chr, hit.strand, (int) hit.pos, (int) hit.pos
                    + hit.matchLength, hit.getK(), cigar, query.toString(), qual, 1, numHits, null);
            return rec;
        }
        else {
            ReadHit split = hit.nextSplit;
            ACGTSequence s1 = query.subSequence(0, hit.matchLength);
            ACGTSequence s2 = query.subSequence(hit.matchLength, m);
            int b1 = Math.min(qual.length(), hit.matchLength);
            int b2 = Math.min(qual.length(), m);
            String q1 = qual.substring(0, b1);
            String q2 = qual.substring(b1, b2);

            if (hit.isUnique()) {
                if (split.isUnique()) {
                    if (hit.chr.equals(split.chr)) {
                        // TODO cigar, score, quality value trimming
                        CIGAR cigar1 = new CIGAR(hit.cigar);
                        //cigar1.add(hit.matchLength, Type.Matches);
                        // cigar1.add(split.matchLength, Type.HardClip);
                        CIGAR cigar2 = new CIGAR(split.cigar);
                        //cigar2.add(hit.matchLength, Type.HardClip);
                        //cigar2.add(split.matchLength, Type.Matches);
                        AlignmentRecord rec = new AlignmentRecord(read.name(), hit.chr, hit.strand, (int) hit.pos,
                                (int) hit.pos + hit.matchLength, hit.diff, cigar1, s1.toString(), q1, 1, 1, null);
                        AlignmentRecord splitRec = new AlignmentRecord(read.name(), split.chr, split.strand,
                                (int) split.pos, (int) split.pos + split.matchLength, split.diff, cigar2,
                                s2.toString(), q2, 1, 1, null);
                        rec.split = splitRec;
                        return rec;
                    }
                    else {
                        // use longer alignment as a base
                        if (hit.matchLength >= split.matchLength) {
                            CIGAR cigar = new CIGAR(hit.cigar);
                            //cigar.add(hit.matchLength, Type.Matches);
                            cigar.add(split.matchLength, Type.SoftClip);
                            AlignmentRecord rec = new AlignmentRecord(read.name(), hit.chr, hit.strand, (int) hit.pos,
                                    (int) hit.pos + m, hit.diff, cigar, query.toString(), qual, 1, numHits, null);
                            return rec;
                        }
                        else {
                            // TODO Use soft clip to represent split part
                            CIGAR cigar = new CIGAR(String.format("%dS%s", hit.matchLength, split.cigar));
                            //cigar.add(hit.matchLength, Type.SoftClip);
                            //cigar.add(split.matchLength, Type.Matches);
                            AlignmentRecord rec = new AlignmentRecord(read.name(), split.chr, split.strand,
                                    (int) split.pos - hit.matchLength, (int) split.pos - hit.matchLength + m,
                                    split.diff, cigar, query.toString(), qual, 1, numHits, null);
                            return rec;
                        }
                    }
                }
                else {
                    CIGAR cigar = new CIGAR(String.format("%s%dS", hit.cigar, split.matchLength));
                    //cigar.add(hit.matchLength, Type.Matches);
                    //cigar.add(split.matchLength, Type.SoftClip);
                    AlignmentRecord rec = new AlignmentRecord(read.name(), hit.chr, hit.strand, (int) hit.pos,
                            (int) hit.pos + m, hit.diff, cigar, query.toString(), qual, 1, numHits, null);
                    return rec;
                }
            }
            else {
                if (split.isUnique()) {
                    // TODO Use soft clip to represent split part
                    CIGAR cigar = new CIGAR(String.format("%dS%s", hit.matchLength, split.cigar));
                    //cigar.add(hit.matchLength, Type.SoftClip);
                    //cigar.add(split.matchLength, Type.Matches);
                    AlignmentRecord rec = new AlignmentRecord(read.name(), split.chr, split.strand, (int) split.pos
                            - hit.matchLength, (int) split.pos - hit.matchLength + m, split.diff, cigar,
                            query.toString(), qual, 1, numHits, null);
                    return rec;
                }

            }

        }
        return null;
    }

}
