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
// IUPACSequenceTest.java
// Since: 2011/04/25
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import static org.junit.Assert.*;

import java.util.Random;

import org.junit.Test;
import org.utgenome.gwt.utgb.client.bio.IUPAC;
import org.xerial.util.log.Logger;

public class IUPACSequenceTest
{
    private static Logger _logger = Logger.getLogger(IUPACSequenceTest.class);

    private final String  orig    = "YGYYYCGCADBACGTNKRWVHT";

    @Test
    public void constructor() throws Exception {
        IUPACSequence s = new IUPACSequence(orig);
        assertEquals(orig.length(), s.textSize());
        for (int i = 0; i < orig.length(); ++i) {
            IUPAC code = s.getIUPAC(i);
            assertEquals(IUPAC.find(String.valueOf(orig.charAt(i))), code);
        }
    }

    @Test
    public void reverse() throws Exception {
        IUPACSequence s = new IUPACSequence(orig);
        IUPACSequence rev = s.reverse();
        for (int i = 0; i < s.textSize(); ++i) {
            IUPAC codeRev = rev.getIUPAC(i);
            IUPAC code = s.getIUPAC(s.textSize() - i - 1);
            assertEquals(code, codeRev);
        }
    }

    @Test
    public void complement() throws Exception {
        IUPACSequence s = new IUPACSequence(orig);
        IUPACSequence r = s.complement();

        for (int i = 0; i < s.textSize(); ++i) {
            IUPAC c = s.getIUPAC(i);
            IUPAC cr = r.getIUPAC(i);
            assertEquals(c.complement(), cr);
        }
    }

    @Test
    public void fastCount() throws Exception {
        Random r = new Random(0);
        StringBuilder seq = new StringBuilder();
        for (int i = 0; i < 69; ++i) {
            seq.append(IUPAC.decode((byte) (r.nextInt(15) + 1)).toString());
        }

        IUPACSequence s = new IUPACSequence(seq.toString());
        for (IUPAC c : IUPAC.values()) {
            for (int x = 0; x < s.textSize(); x++) {
                for (int y = x; y < s.textSize(); y++)
                    assertEquals(String.format("code:%s, s=%d, e=%d", c, x, y), s.count(c, x, y), s.fastCount(c, x, y));
            }
        }

    }

}
