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
// UInt32SAISTest.java
// Since: 2011/02/24
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.sais;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;
import org.xerial.util.log.Logger;

public class UInt32SAISTest
{
    private static Logger _logger = Logger.getLogger(UInt32SAISTest.class);

    @Test
    public void sais() throws Exception {
        LStringSeq s = new LStringSeq("mmiissiissiippii");
        UInt32Array SA = UInt32SAIS.SAIS(s, 65535);

        List<Long> SA_v = new ArrayList<Long>();
        for (long each : SA)
            SA_v.add(each);

        long[] answer = { 16, 15, 14, 10, 6, 2, 11, 7, 3, 1, 0, 13, 12, 9, 5, 8, 4 };
        List<Long> ans = new ArrayList<Long>();
        for (long each : answer)
            ans.add(each);

        _logger.debug(SA_v);
        assertEquals(ans, SA_v);

    }

}
