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
        //BitParallelSmithWaterman.align64(new ACGTSequence(ref), new ACGTSequence("ACGTGGTCTT"));

        //BitParallelSmithWaterman.align64(new ACGTSequence("AAATTT"), new ACGTSequence("AATACTTT"));
        BitParallelSmithWaterman.align64(new ACGTSequence(ref), new ACGTSequence("ACGTGGT"), 2);
        BitParallelSmithWaterman.align64(new ACGTSequence(ref), new ACGTSequence("CTT"), 2);
        //BitParallelSmithWaterman.align64(new ACGTSequence("TATAATAATA"), new ACGTSequence("TAATA"));

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
        //BitParallelSmithWaterman.align64(ref, new ACGTSequence("GGTTCC"), 0);
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
        final int N = 10000;
        StopWatch s1 = new StopWatch();
        s1.stop();
        StopWatch s2 = new StopWatch();
        s2.stop();
        StopWatch s3 = new StopWatch();
        s3.stop();

        for (int k = 0; k < K; ++k) {
            {
                s1.resume();
                for (int i = 0; i < N; ++i) {
                    SmithWatermanAligner.align(ref, query);
                }
                s1.stop();
            }

            {
                s3.resume();
                for (int i = 0; i < N; ++i) {
                    BandedSmithWaterman.align(ref, query);
                }
                s3.stop();
            }

            {
                s2.resume();
                for (int i = 0; i < N; ++i) {
                    BitParallelSmithWaterman.alignBlock(ref, query, (int) query.textSize());
                }
                s2.stop();
            }
        }

        double swTime = s1.getElapsedTime();
        double bpTime = s2.getElapsedTime();
        double bswTime = s3.getElapsedTime();
        _logger.debug("SW: %.2f", swTime);
        _logger.debug("Banded SW: %.2f", bswTime);
        _logger.debug("BP: %.2f", bpTime);
        _logger.debug("SW/BP: %.2f speed up", swTime / bpTime);

    }
}
