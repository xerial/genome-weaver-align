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

import static org.junit.Assert.assertTrue;

import org.junit.Test;
import org.utgenome.weaver.align.BitParallelSmithWaterman.SWResult;
import org.xerial.lens.SilkLens;
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

}
