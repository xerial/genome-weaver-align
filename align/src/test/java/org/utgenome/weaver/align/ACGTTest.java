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
// ACGTTest.java
// Since: 2011/07/13
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import static org.junit.Assert.*;

import org.junit.Test;

public class ACGTTest
{
    @Test
    public void acgt() throws Exception {
        assertEquals(0, ACGT.A.code);
        assertEquals(1, ACGT.C.code);
        assertEquals(2, ACGT.G.code);
        assertEquals(3, ACGT.T.code);
        assertEquals(4, ACGT.N.code);
    }

    @Test
    public void complement() throws Exception {
        assertEquals(ACGT.T, ACGT.A.complement());
        assertEquals(ACGT.G, ACGT.C.complement());
        assertEquals(ACGT.C, ACGT.G.complement());
        assertEquals(ACGT.A, ACGT.T.complement());
        assertEquals(ACGT.N, ACGT.N.complement());

    }
}
