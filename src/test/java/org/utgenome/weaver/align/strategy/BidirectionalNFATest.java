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
// BidirectionalNFATest.java
// Since: 2011/08/08
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.strategy;

import org.junit.BeforeClass;
import org.junit.Test;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.xerial.util.log.Logger;

public class BidirectionalNFATest
{
    private static Logger          _logger = Logger.getLogger(BidirectionalNFATest.class);

    private static FMIndexOnGenome fmIndex;

    @BeforeClass
    public static void setUp() {
        fmIndex = FMIndexOnGenome.buildFromSequence("seq", "AAGCCTAGTTTCCTTG");
    }

    @Test
    public void oneMismatch() throws Exception {
        ACGTSequence q = new ACGTSequence("GCGTAG");
        BidirectionalNFA nfa = new BidirectionalNFA(fmIndex, q);
        nfa.align();
    }

    @Test
    public void forwardExact() throws Exception {
        BidirectionalNFA nfa = new BidirectionalNFA(fmIndex, new ACGTSequence("GCCTA"));
        nfa.align();
    }

    @Test
    public void reverseExact() throws Exception {
        BidirectionalNFA nfa = new BidirectionalNFA(fmIndex, new ACGTSequence("GCCTA").reverseComplement());
        nfa.align();
    }

    @Test
    public void splitExact() throws Exception {
        BidirectionalNFA nfa = new BidirectionalNFA(fmIndex, new ACGTSequence("AAGCCTTCCTTG"));
        nfa.align();
    }

}
