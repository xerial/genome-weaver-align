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
// LLongArrayTest.java
// Since: 2011/02/17
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import org.junit.Test;
import org.xerial.util.StopWatch;
import org.xerial.util.log.Logger;

public class LLongArrayTest
{
    private static Logger _logger = Logger.getLogger(LLongArrayTest.class);

    @Test
    public void performance() throws Exception {

        final int N = 10000000;
        final int K = 10;
        LLongArray ll = new LLongArray(N);
        long[] raw = new long[N];
        StopWatch timer = new StopWatch();

        for (int k = 0; k < K; ++k)
            for (int i = 0; i < N; ++i) {
                raw[i] = i;
            }
        _logger.info("long: " + timer.getElapsedTime());
        timer.reset();
        for (int k = 0; k < K; ++k)
            for (int i = 0; i < N; ++i) {
                ll.set(i, i);
            }
        _logger.info("LLongArray: " + timer.getElapsedTime());

    }
}
