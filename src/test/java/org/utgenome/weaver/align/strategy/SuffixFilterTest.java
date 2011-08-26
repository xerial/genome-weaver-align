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

    private static ACGTSequence         ref     = new ACGTSequence("AAGCCTAGTTTCCTTG");

    @BeforeClass
    public static void setUp() {
        fmIndex = FMIndexOnGenome.buildFromSequence("seq", ref);
    }

    public static void align(String query) throws Exception {
        align(new ACGTSequence(query));
    }

    public static void align(ACGTSequence q) throws Exception {
        SuffixFilter f = new SuffixFilter(fmIndex, ref, config);
        f.align(q);
    }

    @Test
    public void oneMismatchAtTail() throws Exception {
        align("GCCTAA");
    }

    @Test
    public void oneMismatch() throws Exception {
        align("GCGTAGTT");
    }

    @Test
    public void oneMismatchReverse() throws Exception {
        align(new ACGTSequence("GCGTAGTT").reverseComplement());
    }

    @Test
    public void twoMismatch() throws Exception {
        align("TTGCCTAGTTT");
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
    public void longRead() throws Exception {
        align("AAGCCTAGATTCCGTG");
    }

    @Test
    public void clip() throws Exception {
        align("AAGCCTAGGGTCTTT");
    }

}
