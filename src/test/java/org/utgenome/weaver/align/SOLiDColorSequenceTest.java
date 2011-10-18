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
// SOLiDColorSequenceTest.java
// Since: 2011/10/18
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import static org.junit.Assert.*;

import org.junit.Test;
import org.xerial.util.log.Logger;

public class SOLiDColorSequenceTest
{
    private static Logger _logger = Logger.getLogger(SOLiDColorSequenceTest.class);

    @Test
    public void constructor() throws Exception {
        String cStr = "T02012121320120101321310213101330030213020200121311";
        SOLiDColorSequence c = new SOLiDColorSequence(cStr);
        _logger.debug(c.toString());
        assertEquals(cStr, c.toString());
        ACGTSequence acgt = c.toACGTSequence();
        _logger.debug("%s length:%d", acgt, acgt.length());

        SOLiDColorSequence c2 = new SOLiDColorSequence(acgt);
        _logger.debug(c2);
        assertEquals(acgt, c2.toACGTSequence());
    }

    @Test
    public void constructor2() throws Exception {
        SOLiDColorSequence c = new SOLiDColorSequence("T30203212333201220322311002131203003300300311022000");
        _logger.debug(c);
        _logger.debug(c.toACGTSequence());
        _logger.debug(c.reverseComplement());
        _logger.debug(c.reverseComplement().toACGTSequence());
    }

    @Test
    public void reverseComplement() throws Exception {
        String cStr = "T02012121320120101321310213101330030213020200121311";
        SOLiDColorSequence c = new SOLiDColorSequence(cStr);
        SOLiDColorSequence r = c.reverseComplement();

        _logger.debug(c.toACGTSequence());
        _logger.debug(c);
        _logger.debug(c.toACGTSequence().reverseComplement());
        _logger.debug(r);
        _logger.debug(r.toACGTSequence());
    }

    @Test
    public void read() throws Exception {
        SOLiDColorSequence c = new SOLiDColorSequence("T00320010320010320010320010320010320010320010320010");
        _logger.debug(c);
        _logger.debug(c.toACGTSequence());
        _logger.debug(c.reverseComplement());
        _logger.debug(c.reverseComplement().toACGTSequence());
        _logger.debug(c.toACGTSequence().reverseComplement());

    }

}
