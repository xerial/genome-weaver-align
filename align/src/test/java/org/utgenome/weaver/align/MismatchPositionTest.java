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
// MismatchPositionTest.java
// Since: 2011/02/17
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import static org.junit.Assert.*;

import org.junit.Test;

public class MismatchPositionTest
{
    @Test
    public void mismatch() throws Exception {
        MismatchPosition m = new MismatchPosition(10);
        MismatchPosition m2 = m.copyAndSet(3);
        m2 = m2.copyAndSet(7);
        m2 = m2.copyAndSet(9);
        assertEquals("0001000101", m2.toString());
        assertEquals("0000000000", m.toString());

        m2 = m2.copyAndReset(3);
        assertEquals("0000000101", m2.toString());
    }
}
