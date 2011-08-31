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
// QueryMaskTest.java
// Since: 2011/08/31
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import static org.junit.Assert.*;

import org.junit.Test;
import org.xerial.util.log.Logger;

public class QueryMaskTest
{
    private static Logger _logger = Logger.getLogger(QueryMaskTest.class);

    public void check(ACGTSequence query, String answer, ACGT ch, int offset, int boundary) {
        int m = query.length();
        QueryMask qm = new QueryMask(query);
        _logger.debug("ans:%s, offset:%d, boundary:%d", answer, offset, boundary);
        assertEquals(BitVector.parseString(answer),
                BitVector.parseLong(qm.getPatternMaskIn64bitForBidirectionalSearch(ch, offset, boundary), m));
    }

    @Test
    public void mask64() throws Exception {

        ACGTSequence query = new ACGTSequence("AAGATTGC");
        // AAG|ATTGC
        //     ATTGC|GAA
        //      TTGC|GAA
        //       TGC|GAA
        check(query, "00001011", ACGT.A, 0, 0);
        check(query, "00000101", ACGT.A, 1, 0);
        check(query, "11000001", ACGT.A, 3, 3);
        check(query, "01100000", ACGT.A, 4, 3);
        check(query, "00110000", ACGT.A, 5, 3);
    }

    @Test
    public void backwardMask() throws Exception {
        ACGTSequence query = new ACGTSequence("GCCAAGTT");
        //        check(query, "01100000", ACGT.C, 4, 4);
        //        check(query, "00110000", ACGT.C, 5, 4);
        //        check(query, "00011000", ACGT.C, 6, 4);
        //        check(query, "00001100", ACGT.C, 7, 4);
        //        check(query, "00000110", ACGT.C, 8, 4);
        check(query, "00000011", ACGT.C, 9, 4);
        check(query, "00000001", ACGT.C, 10, 4);
    }

}
