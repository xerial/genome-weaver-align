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
// KMPMatcherTest.java
// Since: 2011/10/12
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import static org.junit.Assert.*;

import java.util.List;

import org.junit.Test;
import org.xerial.util.log.Logger;

public class KMPMatcherTest
{
    private static Logger _logger = Logger.getLogger(KMPMatcherTest.class);

    @Test
    public void find() throws Exception {

        ACGTSequence pattern = new ACGTSequence("CCGGG");
        ACGTSequence seq = new ACGTSequence("TTGTCCGGGCCGGGATCCGGGT");

        List<Integer> pos = KMPMatcher.find(seq, pattern);
        assertEquals(3, pos.size());
        _logger.debug(pos);
        assertTrue(pos.contains(4));
        assertTrue(pos.contains(9));
        assertTrue(pos.contains(16));
    }
}
