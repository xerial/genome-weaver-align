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

public class SuffixFilterTest
{
    private static Logger          _logger = Logger.getLogger(SuffixFilterTest.class);

    private static FMIndexOnGenome fmIndex;
    private static AlignmentConfig config  = new AlignmentConfig();

    private static ACGTSequence    ref     = new ACGTSequence("AAGCCTAGTTTCCTTG");

    @BeforeClass
    public static void setUp() {
        fmIndex = FMIndexOnGenome.buildFromSequence("seq", ref);
    }

    public static AlignmentRecord align(String query) throws Exception {
        return align(new ACGTSequence(query));
    }

    public static AlignmentRecord align(ACGTSequence q) throws Exception {
        SuffixFilter f = new SuffixFilter(fmIndex, ref, config);
        List<AlignmentRecord> result = f.align(q);
        if (result.size() == 0)
            return null;
        else
            return result.get(0);
    }

    @Test
    public void oneMismatchAtTail() throws Exception {
        AlignmentRecord a = align("GCCTAC");
        assertEquals("5M1S", a.cigar.toCIGARString());
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
        assertEquals("8M", a.cigar.toCIGARString());
        assertEquals(3, a.start);
        assertEquals(Strand.FORWARD, a.strand);
        assertEquals(1, a.numMismatches);
    }

    @Test
    public void oneMismatchReverse() throws Exception {
        AlignmentRecord a = align(new ACGTSequence("GCGTAGTT").reverseComplement());
        assertEquals("8M", a.cigar.toCIGARString());
        assertEquals(3, a.start);
        assertEquals(Strand.REVERSE, a.strand);
        assertEquals(1, a.numMismatches);
    }

    @Test
    public void bidirectionalSearch() throws Exception {
        AlignmentRecord a = align("AACCCTAGTTTCGTT");
        assertEquals("15M", a.cigar.toCIGARString());
        assertEquals(1, a.start);
        assertEquals(Strand.FORWARD, a.strand);
        assertEquals(2, a.numMismatches);
    }

    @Test
    public void twoMismatchAtHead() throws Exception {
        AlignmentRecord a = align("TTGCCTAGTTT");
        assertEquals("2S9M", a.cigar.toCIGARString());
        assertEquals(3, a.start);
        assertEquals(Strand.FORWARD, a.strand);
        assertEquals(0, a.numMismatches);
    }

    @Test
    public void forwardExact() throws Exception {
        align("GCCTAGT");
    }

    @Test
    public void reverseExact() throws Exception {
        align(new ACGTSequence("GCCTAGT").reverseComplement());
    }

    @Test
    public void splitExact() throws Exception {
        align("AAGCCTTCCTTG");
    }

    @Test
    public void splitExactReverse() throws Exception {
        align(new ACGTSequence("AAGCCTTCCTTG").reverseComplement());
    }

    @Test
    public void split2() throws Exception {
        align("AGCCTATTCCTT");
    }

    @Test
    public void longRead() throws Exception {
        align("AAGCCTAGATTCCGTG");
    }

    @Test
    public void clip() throws Exception {
        align("AAGCCTAGGGTCTTT"); // 7S
    }

    @Test
    public void oneDeletion() throws Exception {
        //     AAGCCTAGTTT 
        //     AAGCCT-GTTT  6M1D4M
        align("AAGCCTGTTT");
    }

    @Test
    public void oneInsertion() throws Exception {
        //     AAGCCT-AGTTT 
        //     AAGCCTTAGTTT  6M1I5M
        align("AAGCCTTAGTTT");
    }

    @Test
    public void twoInsertion() throws Exception {
        //     AAGCCT--AGTTT 
        //     AAGCCTGGAGTTT  6M2I5M
        align("AAGCCTGGAGTTT");
    }

}
