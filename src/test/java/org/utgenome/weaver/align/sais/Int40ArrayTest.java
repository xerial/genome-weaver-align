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
// Int40ArrayTest.java
// Since: 2011/06/29
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

public class Int40ArrayTest
{
    private static Logger _logger = Logger.getLogger(Int40ArrayTest.class);

    @Test
    public void rangeCheck() throws Exception {
        final long max = Int40Array.MAX_VALUE;

        Int40Array array = new Int40Array(1);
        for (long i = max - 10; i <= max; ++i) {
            array.set(0, i);
            assertEquals(i, array.lookup(0));
        }

        array.set(0, -1);
        assertEquals(-1, array.lookup(0));
        array.set(0, -2);
        assertEquals(-2, array.lookup(0));

        final long min = Int40Array.MIN_VALUE;
        for (long i = min + 10; i >= min; --i) {
            array.set(0, i);
            assertEquals(i, array.lookup(0));
        }

    }

    @Test
    public void exhaustiveTest() throws Exception {
        int N = 100;
        Int40Array array = new Int40Array(N);

        long split = (Int40Array.MAX_VALUE - Int40Array.MIN_VALUE) / 1001;
        for (long index = 0; index < N; ++index) {
            for (long i = Int40Array.MIN_VALUE; i <= Int40Array.MAX_VALUE; i += split) {
                array.set(index, i);
                assertEquals(String.format("array[%d]", index), i, array.lookup(index));
            }
        }
    }

    @Test
    public void performance() throws Exception {
        StopWatch timer = new StopWatch();
        int N = 10000000;

        {
            timer.reset();
            long[] u = new long[10];
            for (long i = 0; i < N; ++i) {
                u[(int) (i % 10)] = i;
                long val = u[(int) i % 10];
            }
            _logger.info("long primitive array: " + timer.getElapsedTime());
        }

        {
            timer.reset();
            Int40Array u = new Int40Array(10);
            for (long i = 0; i < N; ++i) {
                u.set(i % 10, i);
                u.lookup(i % 10);
            }
            _logger.info("Int40Array: " + timer.getElapsedTime());
        }

        {
            timer.reset();
            LIntArray array = new LIntArray(10);
            for (long i = 0; i < N; ++i) {
                array.set(i % 10, i);
                array.lookup(i % 10);
            }
            _logger.info("LIntArray: " + timer.getElapsedTime());
        }

    }

}
