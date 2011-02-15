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
                sparseSA[i] = readInt(in);
            }
            return new SparseSuffixArray(sparseSA, N, L);
        }
        finally {
            if (in != null)
                in.close();
        }
    }

    public static SparseSuffixArray buildFromSuffixArray(int[] SA, int L) {
        int sparseSA_length = (SA.length + L - 1) / L;
        int[] sparseSA = new int[sparseSA_length];
        for (int i = 0; i < sparseSA_length; ++i) {
            sparseSA[i] = SA[i * L];
        }
        return new SparseSuffixArray(sparseSA, SA.length, L);
    }

    //    public static SparseSuffixArray buildFromFMIndex(FMIndex fmIndex, int L) {
    //        // TODO use long[] instead of int[]
    //        int N = (int) fmIndex.textSize();
    //        int sparseSA_length = ((N + L - 1) / L);
    //        int[] sparseSA = new int[sparseSA_length];
    //        for (int i = 0; i < sparseSA_length; ++i) {
    //            sparseSA[i] = -1;
    //        }
    //        sparseSA[0] = N - 1;
    //        for (int i = L; i < N; i += L) {
    //            computeSA(sparseSA, i, N, L, fmIndex);
    //        }
    //        return new SparseSuffixArray(sparseSA, (int) fmIndex.textSize(), L);
    //    }
    //
    //    public static long computeSA(int[] sparseSA, long cursor, int N, int L, FMIndex fmIndex) {
    //        for (long j = 1; j <= N; j++) {
    //            cursor = fmIndex.suffixLink(cursor);
    //            if (cursor == 0) {
    //                return j - 1;
    //            }
    //            if (cursor % L == 0) {
    //                int pos = (int) (cursor / L);
    //                if (sparseSA[pos] != -1)
    //                    return sparseSA[pos] + j;
    //                else {
    //                    sparseSA[pos] = (int) (computeSA(sparseSA, cursor, N, L, fmIndex));
    //                    return sparseSA[pos] + j;
    //                }
    //            }
    //        }
    //        throw new IllegalStateException(String.format("cannot reach here: cursor=%d)", cursor));
    //    }

    public long get(long index, FMIndex fmIndex) {
        long pos = index / L;
        long offset = index % L;

        if (pos > Integer.MAX_VALUE) {
            throw new ArrayIndexOutOfBoundsException("index: " + index);
        }

        if (offset == 0)
            return sparseSA[(int) pos];

        long cursor = index;
        final long N = fmIndex.textSize();
        for (long j = 1; j <= N; j++) {
            cursor = fmIndex.suffixLink(cursor);
            if (cursor == 0) {
                return j - 1;
            }
            if (cursor % L == 0)
                return sparseSA[(int) cursor / L] + j;
        }
        throw new IllegalStateException(String.format("cannot reach here: get(index:%d)", index));

    }
}
