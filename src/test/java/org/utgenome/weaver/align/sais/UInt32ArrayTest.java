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
// UInt32ArrayTest.java
// Since: 2011/02/17
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.sais;

import static org.junit.Assert.*;

import org.junit.Test;
import org.utgenome.weaver.align.LIntArray;
import org.xerial.util.StopWatch;
import org.xerial.util.log.Logger;

public class UInt32ArrayTest
{
    private static Logger _logger = Logger.getLogger(UInt32ArrayTest.class);

    @Test
    public void rangeCheck() throws Exception {
        final long max = 0xFFFFFFFFL;

        UInt32Array array = new UInt32Array(1);
        for (long i = max - 10; i <= max; ++i) {
            array.set(0, i);
            assertEquals(i, array.lookup(0));
        }

        array.set(0, -1);
        assertEquals(0xFFFFFFFFL, array.lookup(0));
    }

    @Test
    public void performance() throws Exception {
        StopWatch timer = new StopWatch();
        int N = 10000000;

        timer.reset();
        UInt32Array u = new UInt32Array(10);
        for (long i = 0; i < N; ++i) {
            u.set(i % 10, i);
        }
        _logger.info("UInt32Array: " + timer.getElapsedTime());

        timer.reset();
        LIntArray array = new LIntArray(10);
        for (long i = 0; i < N; ++i) {
            array.set(i % 10, i);
        }
        _logger.info("LIntArray: " + timer.getElapsedTime());

    }

}
