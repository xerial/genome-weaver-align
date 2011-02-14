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
// SparseSuffixArray.java
// Since: 2011/02/12
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

/**
 * 
 * 
 * S(i) = S(X^{(j)}(i)) + j
 * 
 * X(i) = C(BWT[i]) + O(BWT[i], i)
 * 
 * @author leo
 * 
 */
public class SparseSuffixArray
{
    private final int[] sparseSA;
    private final int   N;
    private final int   L;

    private SparseSuffixArray(int[] sparseSA, int N, int L) {
        this.sparseSA = sparseSA;
        this.N = N;
        this.L = L;
    }

    private static void writeInt(int value, OutputStream out) throws IOException {
        out.write((value >>> 24) & 0xFF);
        out.write((value >>> 16) & 0xFF);
        out.write((value >>> 8) & 0xFF);
        out.write((value & 0xFF));
    }

    private static int readInt(InputStream in) throws IOException {
        int value = in.read();
        value <<= 8;
        value |= in.read();
        value <<= 8;
        value |= in.read();
        value <<= 8;
        value |= in.read();
        return value;
    }

    public void saveTo(OutputStream out) throws IOException {
        writeInt(N, out);
        writeInt(L, out);
        writeInt(sparseSA.length, out);
        for (int i = 0; i < sparseSA.length; ++i) {
            writeInt(sparseSA[i], out);
        }
        out.flush();
    }

    public static SparseSuffixArray loadFrom(InputStream in) throws IOException {
        try {
            final int N = readInt(in);
            final int L = readInt(in);
            final int sparseSALength = readInt(in);
            int sparseSA[] = new int[sparseSALength];
            for (int i = 0; i < sparseSA.length; ++i) {
                sparseSA[0] = readInt(in);
            }
            return new SparseSuffixArray(sparseSA, N, L);
        }
        finally {
            if (in != null)
                in.close();
        }
    }

    public static SparseSuffixArray buildFromSuffixArray(int[] SA, int L) {
        int sparseSA_length = SA.length / L + 1;
        int[] sparseSA = new int[sparseSA_length];
        for (int i = 0; i < sparseSA_length; ++i) {
            sparseSA[i] = SA[i * L];
        }
        return new SparseSuffixArray(sparseSA, SA.length, L);
    }

    public int get(int index, FMIndex fmIndex) {
        int pos = index / L;
        int offset = index % L;
        if (offset == 0)
            return sparseSA[pos];

        int cursor = index;
        final int N = fmIndex.textSize();
        for (int j = 1; j <= Integer.MAX_VALUE; j++) {
            cursor = fmIndex.inverseSA(cursor);
            if (cursor == 0) {
                return j - 1;
            }
            if (cursor % L == 0)
                return sparseSA[cursor / L] + j;
        }
        throw new IllegalStateException("cannot reach here");

    }
}
