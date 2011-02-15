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
// BitVectorTest.java
// Since: 2011/02/15
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import static org.junit.Assert.*;

import java.util.TreeSet;

import org.junit.Test;
import org.xerial.util.log.Logger;

public class BitVectorTest
{
    private static Logger _logger = Logger.getLogger(BitVectorTest.class);

    @Test
    public void vector() throws Exception {
        BitVector v = new BitVector(100);
        assertEquals(100, v.size());
        TreeSet<Long> pos = new TreeSet<Long>();
        pos.add(3L);
        pos.add(63L);
        pos.add(80L);
        pos.add(10L);
        pos.add(34L);
        for (long each : pos) {
            v.set(each);
        }
        _logger.debug(v);

        for (long i = 0; i < v.size(); ++i) {
            if (pos.contains(i))
                assertTrue(v.get(i));
            else
                assertFalse(v.get(i));
        }

        long count = 0;
        for (long i = 0; i < v.size(); ++i) {
            assertEquals("index i=" + i, count, v.rank(true, i));
            assertEquals("index i=" + i, i - count, v.rank(false, i));
            if (pos.contains(i))
                count++;
        }
    }

    @Test
    public void select() throws Exception {
        BitVector v = new BitVector(100);
        assertEquals(100, v.size());
        TreeSet<Long> pos = new TreeSet<Long>();
        pos.add(3L);
        pos.add(10L);
        pos.add(34L);
        pos.add(63L);
        pos.add(80L);
        for (long each : pos) {
            v.set(each);
        }
        _logger.debug(v);

        long count = 0;
        for (long each : pos) {
            count++;
            assertEquals(each, v.select(true, count));
        }

    }
}
