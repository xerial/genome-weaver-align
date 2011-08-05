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

import org.junit.Test;
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
        BitParallelSmithWaterman.alignBlock(ref, new ACGTSequence("ACGTGGT"), 2);
    }

    @Test
    public void alignBlockOver64() throws Exception {
        ACGTSequence ref = new ACGTSequence(
                "ACGTGGTCTTACGTTATATAGGGCGCGGGGATATCTGCTATAAAATTAAGCGCAGATGGCGAGGCGTTTTAGAGAGCCGAAGAGACGCGCGCTTCGCGAGAAAGCGGAGCGCAGAGAGGGGCGCCATATATATAGAGCGCGAAGAGAGAT");
        BitParallelSmithWaterman.alignBlock(ref, ref, 2);

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
