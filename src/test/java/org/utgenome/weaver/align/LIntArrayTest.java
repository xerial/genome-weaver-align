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
// LIntArrayTest.java
// Since: 2011/02/16
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import static org.junit.Assert.*;

import org.junit.Test;

public class LIntArrayTest
{
    @Test
    public void overflow() throws Exception {
        LIntArray array = new LIntArray(5);
        final long ans = Integer.MAX_VALUE + 10L;
        array.set(0, ans);
        long v = array.get(0);
        assertEquals(ans, v);
    }

    @Test
    public void exhaustiveCheck() throws Exception {
        LIntArray array = new LIntArray(1);
        for (long i = -0xFFFFFFFF; i < 0xFFFFFFFF; ++i) {
            array.set(0, i);
            assertEquals(i, array.get(i));
        }
    }

    @Test
    public void neg() throws Exception {
        LIntArray array = new LIntArray(1);
        long v = ~11L;
        array.set(0, v);
        assertEquals(v, array.get(0));

    }
}
