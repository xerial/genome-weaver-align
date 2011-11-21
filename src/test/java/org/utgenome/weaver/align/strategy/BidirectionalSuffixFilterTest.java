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
// SuffixFilterTest.java
// Since: 2011/07/27
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.strategy;

import static org.junit.Assert.*;

import java.util.List;

import org.junit.BeforeClass;
import org.junit.Test;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.AlignmentConfig;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.record.AlignmentRecord;
import org.xerial.util.log.Logger;

public class BidirectionalSuffixFilterTest
{
    private static Logger          _logger = Logger.getLogger(BidirectionalSuffixFilterTest.class);

    private static FMIndexOnGenome fmIndex;
    private static AlignmentConfig config  = new AlignmentConfig();

    private static ACGTSequence    ref     = new ACGTSequence("AAGCCTAGTTTCCTTG");

    @BeforeClass
    public static void setUp() {
        fmIndex = FMIndexOnGenome.buildFromSequence("seq", ref);
        config.k = 2;
    }

    public static AlignmentRecord align(String query) throws Exception {
        return align(new ACGTSequence(query));
    }

    public static AlignmentRecord align(ACGTSequence q) throws Exception {
        BidirectionalSuffixFilter f = new BidirectionalSuffixFilter(fmIndex, ref, config);
        List<AlignmentRecord> result = f.align(q);
        if (result.size() == 0)
            return null;
        else
            return result.get(0);
    }

    @Test
    public void oneMismatchAtTail() throws Exception {
        AlignmentRecord a = align("GCCTAC");
        assertEquals("5M1S", a.cigar.toString());
        assertEquals(3, a.start);
        assertEquals(Strand.FORWARD, a.strand);
        assertEquals(0, a.numMismatches);
    }

    @Test
    public void oneMismatch() throws Exception {
        // GCCTAGTT
        // |||X||||
        // GCCAAGTT
        AlignmentRecord a = align("GCCAAGTT");
        assertEquals("8M", a.cigar.toString());
        assertEquals(3, a.start);
        assertEquals(Strand.FORWARD, a.strand);
        assertEquals(1, a.numMismatches);
    }

    @Test
    public void oneMismatchReverse() throws Exception {
        AlignmentRecord a = align(new ACGTSequence("GCGTAGTT").reverseComplement());
        assertEquals("8M", a.cigar.toString());
        assertEquals(3, a.start);
        assertEquals(Strand.REVERSE, a.strand);
        assertEquals(1, a.numMismatches);
    }

    @Test
    public void bidirectionalSearch() throws Exception {
        AlignmentRecord a = align("AACCCTAGTTTCGTT");
        assertEquals("15M", a.cigar.toString());
        assertEquals(1, a.start);
        assertEquals(Strand.FORWARD, a.strand);
        assertEquals(2, a.numMismatches);
    }

    @Test
    public void bidirectionalSearchReverse() throws Exception {
        AlignmentRecord a = align(new ACGTSequence("AACCCTAGTTTCGTT").reverseComplement());
        assertEquals("15M", a.cigar.toString());
        assertEquals(1, a.start);
        assertEquals(Strand.REVERSE, a.strand);
        assertEquals(2, a.numMismatches);
    }

    @Test
    public void twoMismatchAtHead() throws Exception {
        AlignmentRecord a = align("TTGCCTAGTTT");
        assertEquals("2S9M", a.cigar.toString());
        assertEquals(3, a.start);
        assertEquals(Strand.FORWARD, a.strand);
        assertEquals(0, a.numMismatches);
    }

    @Test
    public void forwardExact() throws Exception {
        AlignmentRecord a = align("GCCTAGT");
        assertEquals("7M", a.cigar.toString());
        assertEquals(3, a.start);
        assertEquals(Strand.FORWARD, a.strand);
        assertEquals(0, a.numMismatches);
    }

    @Test
    public void reverseExact() throws Exception {
        AlignmentRecord a = align(new ACGTSequence("GCCTAGT").reverseComplement());
        assertEquals("7M", a.cigar.toString());
        assertEquals(3, a.start);
        assertEquals(Strand.REVERSE, a.strand);
        assertEquals(0, a.numMismatches);
    }

    @Test
    public void splitExact() throws Exception {
        AlignmentRecord a = align("AAGCCTATCCTTG");
        assertEquals("7M", a.cigar.toString());
        assertEquals(1, a.start);
        assertEquals(Strand.FORWARD, a.strand);
        assertEquals(0, a.numMismatches);
        AlignmentRecord s = a.split;
        assertNotNull(s);
        assertEquals("6M", s.cigar.toString());
        assertEquals(11, s.start);
        assertEquals(Strand.FORWARD, s.strand);
        assertEquals(0, s.numMismatches);

    }

    @Test
    public void splitExactReverse() throws Exception {
        AlignmentRecord a = align(new ACGTSequence("AAGCCTATCCTTG").reverseComplement());
        assertEquals("7M", a.cigar.toString());
        assertEquals(1, a.start);
        assertEquals(Strand.REVERSE, a.strand);
        assertEquals(0, a.numMismatches);
        AlignmentRecord s = a.split;
        assertNotNull(s);
        assertEquals("6M", s.cigar.toString());
        assertEquals(11, s.start);
        assertEquals(Strand.REVERSE, s.strand);
        assertEquals(0, s.numMismatches);

    }

    @Test
    public void split2() throws Exception {
        align("AGCCTATTCCTT");
    }

    @Test
    public void longRead() throws Exception {
        AlignmentRecord a = align("AAGCCTAGATTCCGTG");
        assertEquals(1, a.start);
        assertEquals(Strand.FORWARD, a.strand);
        assertEquals(2, a.numMismatches);
        assertEquals("16M", a.cigar.toString());
    }

    @Test
    public void clip() throws Exception {
        // r:AAGCCTAGTTTCCTTG
        //   ||||||||||           
        // q:AAGCCTAGTTAAAAAA
        AlignmentRecord a = align("AAGCCTAGTTAAAAAA");
        assertEquals(1, a.start);
        assertEquals(11, a.end);
        assertEquals("10M6S", a.cigar.toString());
    }

    @Test
    public void clip2() throws Exception {
        // r: AAGCCTAGTTTCCTTG
        //          ||||||||||
        // q:TTTTTTGAGTTTCCTTG
        AlignmentRecord a = align("TTTTTTGAGTTTCCTTG");
        assertEquals(7, a.start);
        assertEquals(16, a.end);
        assertEquals("7S10M", a.cigar.toString());
    }

    @Test
    public void oneDeletion() throws Exception {
        //     AAGCCTAGTTT 
        //     AAGCCT-GTTT  6M1D4M
        AlignmentRecord a = align("AAGCCTGTTT");
        assertEquals(1, a.start);
        assertEquals(Strand.FORWARD, a.strand);
        assertEquals(1, a.numMismatches);
        assertEquals("6M1D4M", a.cigar.toString());
    }

    @Test
    public void oneInsertion() throws Exception {
        //     AAGCCT-AGTT
        //     AAGCCTCAGTT  6M1I4M
        AlignmentRecord a = align("AAGCCTCAGTT");
        assertEquals(1, a.start);
        assertEquals(Strand.FORWARD, a.strand);
        assertEquals(1, a.numMismatches);
        assertEquals("6M1I4M", a.cigar.toString());
    }

    @Test
    public void twoInsertion() throws Exception {
        //     AAGCC--TAGTTT
        //     AAGCCAATAGTTT  5M2I7M
        AlignmentRecord a = align("AAGCCAATAGTTT");
        assertEquals(1, a.start);
        assertEquals(Strand.FORWARD, a.strand);
        assertEquals(2, a.numMismatches);
        assertEquals("5M2I6M", a.cigar.toString());
    }

}
