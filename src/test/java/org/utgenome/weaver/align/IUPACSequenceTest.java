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

import org.junit.Test;
import org.utgenome.gwt.utgb.client.bio.IUPAC;

public class IUPACSequenceTest
{
    private final String orig = "YGYYYCGCADBACGTNKRWVHT";

    @Test
    public void constructor() throws Exception {
        IUPACSequence s = new IUPACSequence(orig);
        assertEquals(orig.length() + 1, s.textSize());
        for (int i = 0; i < orig.length(); ++i) {
            IUPAC code = s.getIUPAC(i);
            assertEquals(IUPAC.find(String.valueOf(orig.charAt(i))), code);
        }
    }

    @Test
    public void reverse() throws Exception {
        IUPACSequence s = new IUPACSequence(orig);
        IUPACSequence rev = s.reverse();
        for (int i = 0; i < s.textSize() - 1; ++i) {
            IUPAC codeRev = rev.getIUPAC(i);
            IUPAC code = s.getIUPAC(s.textSize() - i - 2);
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

}
