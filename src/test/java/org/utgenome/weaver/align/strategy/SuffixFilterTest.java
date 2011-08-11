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

import org.junit.BeforeClass;
import org.junit.Test;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.AlignmentScoreConfig;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.xerial.util.log.Logger;

public class SuffixFilterTest
{
    private static Logger               _logger = Logger.getLogger(SuffixFilterTest.class);

    private static FMIndexOnGenome      fmIndex;
    private static AlignmentScoreConfig config  = new AlignmentScoreConfig();

    @BeforeClass
    public static void setUp() {
        fmIndex = FMIndexOnGenome.buildFromSequence("seq", "AAGCCTAGTTTCCTTG");
    }

    @Test
    public void oneMismatch() throws Exception {
        ACGTSequence q = new ACGTSequence("GCGTAG");
        SuffixFilter f = new SuffixFilter(fmIndex, config, q.textSize());
        f.align(q);
    }

    @Test
    public void twoMismatch() throws Exception {
        ACGTSequence q = new ACGTSequence("TTGCCTAGTTT");
        SuffixFilter f = new SuffixFilter(fmIndex, config, q.textSize());
        f.align(q);
    }

    @Test
    public void forwardExact() throws Exception {
        ACGTSequence q = new ACGTSequence("GCCTAGT");
        SuffixFilter f = new SuffixFilter(fmIndex, config, q.textSize());
        f.align(q);
    }

    @Test
    public void reverseExact() throws Exception {
        ACGTSequence q = new ACGTSequence("GCCTAGT").reverseComplement();
        SuffixFilter f = new SuffixFilter(fmIndex, config, q.textSize());
        f.align(q);
    }

    @Test
    public void splitExact() throws Exception {
        ACGTSequence q = new ACGTSequence("AAGCCTTCCTTG");
        SuffixFilter f = new SuffixFilter(fmIndex, config, q.textSize());
        f.align(q);
    }

}
