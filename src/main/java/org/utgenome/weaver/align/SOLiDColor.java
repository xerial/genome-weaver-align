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
// SOLiDColor.java
// Since: 2011/07/20
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.util.ArrayList;
import java.util.List;

import org.xerial.util.CollectionUtil;
import org.xerial.util.Functor;
import org.xerial.util.StringUtil;

/**
 * Color code sequence character for SOLiD
 * 
 * @author leo
 * 
 */
public enum SOLiDColor {
    C0(0), C1(1), C2(2), C3(3), N(4);

    public final int                  code;

    private final static SOLiDColor[] codeTable            = { C0, C1, C2, C3, N, N, N, N, N };
    private final static String[]     codeName             = { "0", "1", "2", "3", "." };

    private final static SOLiDColor[] charToCode           = { N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N,
            N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, C0, C1, C2, C3, N,
            N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N,
            N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N,
            N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N,
            N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N,
            N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N,
            N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N };

    private final static ACGT[]       base_color2baseTable = {
                                                           // A0 A1 A2 A3

                                                           };

    private SOLiDColor(int code) {
        this.code = code;
    }

    @Override
    public String toString() {
        return codeName[code];
    }

    public static SOLiDColor decode(int value) {
        return codeTable[value & 0x07];
    }

    public static SOLiDColor encode(char ch) {
        return charToCode[ch & 0xFF];
    }

    public static List<SOLiDColor> genCharToColorCodeTable() {

        ArrayList<SOLiDColor> table = new ArrayList<SOLiDColor>();
        for (int i = 0; i < 255; ++i) {
            char c = (char) i;
            SOLiDColor code = SOLiDColor.N;
            switch (c) {
            case '0':
                code = SOLiDColor.C0;
                break;
            case '1':
                code = SOLiDColor.C1;
                break;
            case '2':
                code = SOLiDColor.C2;
                break;
            case '3':
                code = SOLiDColor.C3;
                break;
            }
            table.add(code);
        }

        return table;
    }
}
