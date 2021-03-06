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
// BitParallelSmithWatermanTest.java
// Since: 2011/07/29
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import static org.junit.Assert.*;

import org.junit.Test;
import org.utgenome.weaver.align.SmithWatermanAligner.Alignment;
import org.utgenome.weaver.align.record.SWResult;
import org.xerial.lens.SilkLens;
import org.xerial.util.StopWatch;
import org.xerial.util.log.Logger;

public class BitParallelSmithWatermanTest
{
    private static Logger _logger = Logger.getLogger(BitParallelSmithWatermanTest.class);

    @Test
    public void align64() throws Exception {

        String ref = "ACGTGGTCTT";
        //BitParallelSmithWaterman.align64(new ACGTSeq(ref), new ACGTSeq("ACGTGGTCTT"));

        //BitParallelSmithWaterman.align64(new ACGTSeq("AAATTT"), new ACGTSeq("AATACTTT"));
        BitParallelSmithWaterman.align64(new ACGTSequence(ref), new ACGTSequence("ACGTGGT"), 2);
        BitParallelSmithWaterman.align64(new ACGTSequence(ref), new ACGTSequence("CTT"), 2);
        //BitParallelSmithWaterman.align64(new ACGTSeq("TATAATAATA"), new ACGTSeq("TAATA"));

    }

    @Test
    public void alignBlock() throws Exception {
        ACGTSequence ref = new ACGTSequence("ACGTGGTCTT");
        SWResult r = BitParallelSmithWaterman.alignBlock(ref, new ACGTSequence("ACGTGGT"), 2);
        _logger.debug(SilkLens.toSilk("result", r));
    }

    @Test
    public void alignBlocks() throws Exception {
        ACGTSequence ref = new ACGTSequence("ACGTCATA");
        SWResult r = BitParallelSmithWaterman.alignBlock(ref, ref, 0, 5);
        _logger.debug(SilkLens.toSilk("result", r));
    }

    @Test
    public void alignBlocks2() throws Exception {
        ACGTSequence ref = new ACGTSequence("ACGTCATAACG");
        SWResult r = BitParallelSmithWaterman.alignBlock(ref, ref, 1);
        _logger.debug(SilkLens.toSilk("result", r));
    }

    @Test
    public void noMatch() throws Exception {
        ACGTSequence ref = new ACGTSequence("ACGTCATAACG");
        //BitParallelSmithWaterman.align64(ref, new ACGTSeq("GGTTCC"), 0);
        SWResult r = BitParallelSmithWaterman.alignBlock(ref, new ACGTSequence("GGTTCC"), 0);
        assertTrue(r.diff > 0);
        _logger.debug(SilkLens.toSilk("result", r));
    }

    @Test
    public void alignClip() throws Exception {
        ACGTSequence ref = new ACGTSequence("ACGTCATAACG");
        SWResult r = BitParallelSmithWaterman.alignBlock(ref, new ACGTSequence("TAACG"), 6);
        _logger.debug(SilkLens.toSilk("result", r));
        SWResult r2 = BitParallelSmithWaterman.alignBlock(ref, new ACGTSequence("ACGTC"), 6);
        _logger.debug(SilkLens.toSilk("result", r2));
    }

    @Test
    public void clip() throws Exception {
        String ref = "ACGTGGTCTT";
        BitParallelSmithWaterman.align64(new ACGTSequence(ref), new ACGTSequence("ACGTGNNNNNGTCTT"), 5);
        BitParallelSmithWaterman.localAlign64(new ACGTSequence(ref), new ACGTSequence("ACGTGNNNNNGTCTT"), 5);
    }

    @Test
    public void softClip() throws Exception {
        String ref = "ACGTTAC";
        BitParallelSmithWaterman.align64(new ACGTSequence(ref), new ACGTSequence("GGTTA"), 1);
        BitParallelSmithWaterman.localAlign64(new ACGTSequence(ref), new ACGTSequence("GGTTA"), 1);

        BitParallelSmithWaterman.localAlign64(new ACGTSequence(ref), new ACGTSequence("ACG"), 1);
        BitParallelSmithWaterman.localAlign64(new ACGTSequence(ref), new ACGTSequence("TAC"), 1);

        BitParallelSmithWaterman.align64(new ACGTSequence(ref), new ACGTSequence("ACG"), 1);
        BitParallelSmithWaterman.align64(new ACGTSequence(ref), new ACGTSequence("TAC"), 1);

    }

    @Test
    public void perfComparison() throws Exception {

        ACGTSequence ref = new ACGTSequence(
                "ACGTGCGGGTTATATTACCGGTTCGCGGCATAGGAAATTGGATAACCGCTAGCATGCATGCCCGATTCAGTGGTGTCACATTTGCCGATCACCCTTCGGTCAGTT");
        ACGTSequence query = new ACGTSequence(
                "TATTACCGGTTCGCGGCATAGGAAATTGGAAAACCGCTAGCATGCATGCCCGATTCAGTGGTGTCACATTTGCCGATC");

        final int K = 5;
        final int N = 1000;
        StopWatch s1 = new StopWatch();
        s1.stop();
        StopWatch s2 = new StopWatch();
        s2.stop();
        StopWatch s3 = new StopWatch();
        s3.stop();
        StopWatch s4 = new StopWatch();
        s4.stop();

        final int W = new AlignmentConfig().bandWidth;

        for (int k = 0; k < K; ++k) {
            {
                s1.resume();
                for (int i = 0; i < N; ++i) {
                    SmithWatermanAligner.standardAlign(ref, query);
                }
                s1.stop();
            }

            {
                s3.resume();
                for (int i = 0; i < N; ++i) {
                    SmithWatermanAligner.bandedAlign(ref, query);
                }
                s3.stop();
            }

            {
                s2.resume();
                for (int i = 0; i < N; ++i) {
                    BitParallelSmithWaterman.alignBlock(ref, query, W);
                }
                s2.stop();
            }

            {
                s4.resume();
                for (int i = 0; i < N; ++i) {
                    BitParallelSmithWaterman.alignBlockDetailed(ref, query, k);
                }
                s4.stop();
            }

        }

        double swTime = s1.getElapsedTime();
        double bpTime = s2.getElapsedTime();
        double bswTime = s3.getElapsedTime();
        double bpDetailedTime = s4.getElapsedTime();
        _logger.debug("SW: %.2f", swTime);
        _logger.debug("Banded SW: %.2f", bswTime);
        _logger.debug("BitParallel: %.2f", bpTime);
        _logger.debug("BitParallel Detailed: %.2f", bpDetailedTime);
        _logger.debug("BitParallel Detailed/BitParallel: %.2f speed up", bpDetailedTime / bpTime);
        _logger.debug("SW/BitParallel: %.2f speed up", swTime / bpTime);
        _logger.debug("Banded SW/BitParallel: %.2f speed up", bswTime / bpTime);

    }

    @Test
    public void deletion() throws Exception {
        Alignment alignment = BitParallelSmithWaterman.alignBlockDetailed("TATACCAAGATCTAGAGATCTGG",
                "TACCAAGATAGAGATCTGG", 31);
        _logger.debug(alignment);
        assertEquals(2, alignment.pos);
        assertEquals("8M2D11M", alignment.cigar.toString());
        assertEquals(2, alignment.numMismatches);
    }

    @Test
    public void mismatches() throws Exception {
        Alignment alignment = BitParallelSmithWaterman.alignBlockDetailed("GATCTA", "GCTATA", 31);
        _logger.debug(alignment);
        assertEquals(0, alignment.pos);
        assertEquals("6M", alignment.cigar.toString());
        assertEquals(2, alignment.numMismatches);
    }

    @Test
    public void insertion() throws Exception {
        Alignment alignment = BitParallelSmithWaterman.alignBlockDetailed("TATACCAAGATCTAGAGATCTGG",
                "TACCAAGATCTCTAGAGATCTGG", 31);
        _logger.debug(alignment);
        assertEquals(2, alignment.pos);
        assertEquals("8M2I13M", alignment.cigar.toString());
        assertEquals(2, alignment.numMismatches);
    }

    @Test
    public void local() throws Exception {
        Alignment alignment = BitParallelSmithWaterman.alignBlockDetailed("TATACCAAGATCTAGAGATCTGG", "ACCAAGATCTAGAG",
                31);
        _logger.debug(alignment);
        assertEquals("14M", alignment.cigar.toString());
        assertEquals(3, alignment.pos);
        assertEquals(0, alignment.numMismatches);
    }

    @Test
    public void softClip2() throws Exception {
        Alignment alignment = BitParallelSmithWaterman.alignBlockDetailed("TATACCAAGATCTAGAGATCTGG",
                "GGCGCACCAAGATCTAGAG", 31);
        _logger.debug(alignment);
        assertEquals("5S14M", alignment.cigar.toString());
        assertEquals(3, alignment.pos);
        assertEquals(0, alignment.numMismatches);
    }

    @Test
    public void tailClip() throws Exception {
        Alignment alignment = BitParallelSmithWaterman.alignBlockDetailed("TATACCAAGATCTAGAGTCTGG",
                "ACCAAGATCTAGAGAAAA", 31);
        _logger.debug(alignment);
        assertEquals("14M4S", alignment.cigar.toString());
        assertEquals(3, alignment.pos);
        assertEquals(0, alignment.numMismatches);
    }

    @Test
    public void test1() throws Exception {
        //        Alignment alignment = BitParallelSmithWaterman.alignBlockDetailed(
        //                "CTAGTGTATCTCGTTATTGCTGCGGTTAAAAAGCTCGTAGTTGGATCTAGGTTACGTGCCGCAGTTCGCAATTTGCGTCAACTGTGGTCGTGACTT",
        //                          "GGGTTATTGCTGCGGTTAAAAAGCTCGTAGTTGGATCTAGGTTACGTGCCGCAGTTCGCAATTTGCGTCAACTGTG", 10);

        ACGTSequence r = new ACGTSequence(
                "CTAGTGTATCTCGTTATTGCTG CGGTTATTGC TGCGGTTAAA AAGCTCGTAG TTGGATCTAG GTTACGTGCC GCAGTTCGCA ATTTGCGTCA ACTGTGGTCG TGACTT");

        ACGTSequence q = new ACGTSequence(
                "GGGTTATTGC TGCGGTTAAA AAGCTCGTAG TTGGATCTAG GTTACGTGCC GCAGTTCGCA ATTTGCGTCA");

        Alignment alignment = BitParallelSmithWaterman.alignBlockDetailed(r, q, 10);
        _logger.debug(alignment);

        //        AlignBlocksDetailed aln = new BitParallelSmithWaterman.AlignBlocksDetailed((int) q.textSize(), 10);
        //        SWResult sw = aln.align(r, q);
        //        _logger.debug("tail pos:%d", sw.tailPos);
        //        Alignment traceback = aln.traceback(r, q, 69);
        //        _logger.debug(traceback);

    }

    @Test
    public void test2() throws Exception {

        ACGTSequence r = new ACGTSequence(
                "CTAGTGTATCTCGTTATTGCTGCGGTTAAAAAGCTCGTAGTTGGATCTAGGTTACGTGCCGCAGTTCGCAATTTGCGTCAACTGTGGTCGTGACTT");
        ACGTSequence q = new ACGTSequence(
                "GGGTTATTGCTGCGGTTAAAAAGCTCGTAGTTGGATCTAGGTTACGTGCCGCAGTTCGCAATTTGCGTCAACTGTG");

        Alignment alignment = BitParallelSmithWaterman.alignBlockDetailed(r, q, 15);
        _logger.debug(alignment);

    }

    @Test
    public void insertion2() throws Exception {
        ACGTSequence r = new ACGTSequence("AAGCCTAGTTTCCTT");
        ACGTSequence q = new ACGTSequence("AAGCCTGGAGTTT");

        Alignment alignment = BitParallelSmithWaterman.alignBlockDetailed(r, q, 15);
        _logger.debug(alignment);

    }

    @Test
    public void noMismatch() throws Exception {
        ACGTSequence r = new ACGTSequence(
                "GCTTCAGTTTCCTGACACTTAAAAAAAAAAGAGTTGCTTATTATTTTAATGAGACTAATGCTTACACTCTGAGTTACTTGTAAGGTGATTGGTTACTTTAATGTTATTATAAGTAATTT");
        ACGTSequence q = new ACGTSequence(
                "CTGACACTTAAAAAAAAAAGAGTTGCTTATTATTTTAATGAGACTAATGCTTACACTCTGAGTTACTTGTAAGGTGATTGGTTACTTTAATGTTATT");

        Alignment alignment = BitParallelSmithWaterman.alignBlockDetailed(r, q, 11);
        _logger.debug(alignment);
        assertEquals("97M", alignment.cigar.toString());
        assertEquals(11, alignment.pos);
        assertEquals(0, alignment.numMismatches);

    }

    @Test
    public void exactMatch() throws Exception {

        ACGTSequence r = new ACGTSequence(
                "GCTTCAGTTTCCTGACACTTAAAAAAAAAAGAGTTGCTTATTATTTTAATGAGACTAATGCTTACACTCTGAGTTACTTGTAAGGTGATTGGTTACTTTAATGTTATTATAAGTAATTTATTGGTTACTTTAATGTTATTATAAGTAAT");

        for (int k = 30; k < r.length(); ++k) {
            for (int i = 0; i < r.length() - k; ++i) {
                ACGTSequence q = r.subSequence(i, i + k);
                Alignment alignment = BitParallelSmithWaterman.alignBlockDetailed(r, q, 2);
                assertEquals(i, alignment.pos);
                assertEquals(String.format("%dM", k), alignment.cigar.toString());
                assertEquals(0, alignment.numMismatches);
            }
        }

    }

}
