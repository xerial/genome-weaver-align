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

/**
 * Color code sequence character for SOLiD.
 * 
 * 
 * Color 0 is for no change, 1 for 00 (A) <=> 01 (C), 10 (G) <=> 11 (T) (flip
 * LSB), 2 for 00 (A) <=> 10 (G), 01 (C) <=> 11 (T) (flip MSB), 3 for complement
 * 00 (A) <=> 11 (T), 01 (C) <=> 10 (G).
 * 
 * 
 * @author leo
 * 
 */
public enum SOLiDColor {
    C0(0), C1(1), C2(2), C3(3), N(4);

    public final int                  code;

    private final static SOLiDColor[] codeTable                = { C0, C1, C2, C3, N, N, N, N, N };
    private final static String[]     codeName                 = { "0", "1", "2", "3", "." };

    private final static SOLiDColor[] charToCode               = { N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N,
            N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, C0, C1, C2,
            C3, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N,
            N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N,
            N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N,
            N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N,
            N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N,
            N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N };

    private final static SOLiDColor[] dibaseToColorTable       = {
                                                               // A C G T
            C0, C1, C2, C3, // A
            C1, C0, C3, C2, // C
            C2, C3, C0, C1, // G
            C3, C2, C1, C0                                    // T
                                                               };

    private final static ACGT[]       baseAndColorToColorTable = {
                                                               // C0, C1, C2, C3
            ACGT.A, ACGT.C, ACGT.G, ACGT.T, // A 
            ACGT.C, ACGT.A, ACGT.T, ACGT.G, // C
            ACGT.G, ACGT.T, ACGT.A, ACGT.C, // G
            ACGT.T, ACGT.G, ACGT.C, ACGT.A, // T
                                                               };

    private SOLiDColor(int code) {
        this.code = code;
    }

    public static ACGT decode(ACGT base, SOLiDColor color) {
        if (base == ACGT.N || color == SOLiDColor.N) {
            return ACGT.N;
        }
        else {
            return baseAndColorToColorTable[(base.code << 2) | color.code];
        }
    }

    public static SOLiDColor encode(ACGT prev, ACGT next) {
        if (prev == ACGT.N || next == ACGT.N)
            return SOLiDColor.N;
        else
            return dibaseToColorTable[(prev.code << 2) | next.code];
    }

    @Override
    public String toString() {
        return codeName[code];
    }

    public static SOLiDColor decode(int code) {
        return codeTable[code & 0x07];
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
