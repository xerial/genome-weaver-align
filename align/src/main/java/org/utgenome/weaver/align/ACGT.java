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
// ACGT.java
// Since: 2011/06/30
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

/**
 * 3-bit encoding of ACGT
 * 
 * @author leo
 * 
 */
public enum ACGT {
    A(0x00, 1), C(0x01, 1 << 1), G(0x02, 1 << 2), T(0x03, 1 << 3), N(0x04, 0x0F);

    private final static byte[] charToACGTCodeTable = { 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 };

    private final static ACGT[] codeTable           = { A, C, G, T, N, N, N, N, N };
    private final static char[] charTable           = { 'A', 'C', 'G', 'T', 'N' };
    public final byte           code;
    public final int            bitFlag;
    public final static ACGT[]  exceptN             = { A, C, G, T };

    private ACGT(int code, int bitFlag) {
        this.code = (byte) code;
        this.bitFlag = bitFlag;

    }

    public ACGT complement() {
        return complement(code);
    }

    public char toChar() {
        return charTable[code];
    }

    public static ACGT complement(byte code) {
        return codeTable[(~code & 0x03) | (code & 0x04)];
    }

    public static ACGT decode(byte code) {
        return codeTable[code & 0x07];
    }

    public static ACGT decode(int code) {
        return decode((byte) code);
    }

    public static ACGT encode(char ch) {
        return decode(charToACGTCodeTable[ch & 0xFF]);
    }

    public static byte to3bitCode(char ch) {
        return charToACGTCodeTable[ch & 0xFF];
    }

    public boolean match(ACGT ch) {
        return (this.bitFlag & ch.bitFlag) != 0;
    }

}
