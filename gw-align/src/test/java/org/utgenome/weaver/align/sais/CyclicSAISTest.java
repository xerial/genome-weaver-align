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

import org.junit.Test;
import org.utgenome.weaver.align.IUPACSequence;
import org.utgenome.weaver.align.LSeq;
import org.xerial.util.log.Logger;

public class CyclicSAISTest
{
    private static Logger _logger = Logger.getLogger(CyclicSAISTest.class);

    public static long[] toLongArray(LSeq seq) {
        long[] a = new long[(int) seq.textSize()];
        for (int i = 0; i < seq.textSize(); ++i)
            a[i] = seq.lookup(i);
        return a;
    }

    @Test
    public void sais() throws Exception {
        LStringSeq s = new LStringSeq("mmiissiissiippii");
        LSeq SA = CyclicSAIS.SAIS(s, 255);

        long[] answer = { 16, 15, 14, 10, 6, 2, 11, 7, 3, 1, 0, 13, 12, 9, 5, 8, 4 };
        assertArrayEquals(answer, toLongArray(SA));
    }

    @Test
    public void saisSeq() {
        //IUPACSequence s = new IUPACSequence("ACGTTA ACGTTA");
        IUPACSequence s = new IUPACSequence("ACT ACTA", true);
        LSeq SA = CyclicSAIS.SAIS(s, 16);
        _logger.debug(SA);

        //long answer[] = { 13, 6, 12, 5, 7, 0, 8, 1, 9, 2, 11, 4, 10, 3 };
        long answer[] = { 9, 4, 0, 7, 5, 1, 6, 2, 8, 3 };

        assertArrayEquals(answer, toLongArray(SA));

    }

    @Test
    public void saisInt() {
        UInt32Array T = new UInt32Array(6);
        int[] t = { 3, 2, 2, 3, 1, 0 };
        for (int i = 0; i < t.length; ++i)
            T.set(i, t[i]);

        LSeq SA = CyclicSAIS.SAIS(T, 4);

        long[] ans = { 5, 4, 1, 2, 3, 0 };
        assertArrayEquals(ans, toLongArray(SA));

    }

    @Test
    public void saisTATA() {
        IUPACSequence s = new IUPACSequence("TATAATAATATAATA", true);
        LSeq SA = CyclicSAIS.SAIS(s, 16);
        _logger.debug(SA);

        long answer[] = { 15, 14, 11, 3, 6, 12, 9, 1, 4, 7, 13, 10, 2, 5, 8, 0 };
        assertArrayEquals(answer, toLongArray(SA));
    }

}
