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

import org.junit.Test;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.parallel.Reporter;
import org.xerial.lens.SilkLens;
import org.xerial.util.log.Logger;

public class BidirectionalNFATest
{
    private static Logger _logger = Logger.getLogger(BidirectionalNFATest.class);

    @Test
    public void sample() throws Exception {

        FMIndexOnGenome fmIndex = FMIndexOnGenome.buildFromSequence("seq", "AACCGTA");
        ACGTSequence q = new ACGTSequence("AAGCGTA");
        BidirectionalNFA nfa = new BidirectionalNFA(fmIndex, q, Strand.FORWARD);

        nfa.bidirectionalAlign(new Reporter() {
            @Override
            public void emit(Object result) throws Exception {
                _logger.debug(SilkLens.toSilk(result));
            }
        });
    }
}
