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

import java.io.DataInputStream;
import java.io.DataOutputStream;
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
    private final LIntArray sparseSA;
    private final long      N;
    private final int       L;

    private SparseSuffixArray(LIntArray sparseSA, long N, int L) {
        this.sparseSA = sparseSA;
        this.N = N;
        this.L = L;
    }

    public void saveTo(OutputStream out) throws IOException {
        DataOutputStream d = new DataOutputStream(out);
        d.writeLong(N);
        d.writeInt(L);
        d.writeLong(sparseSA.size());
        for (int i = 0; i < sparseSA.size(); ++i) {
            d.writeInt((int) sparseSA.get(i));
        }
        d.flush();
    }

    public static SparseSuffixArray loadFrom(InputStream in) throws IOException {
        DataInputStream d = new DataInputStream(in);
        try {
            final long N = d.readLong();
            final int L = d.readInt();
            final long sparseSALength = d.readLong();
            LIntArray sparseSA = new LIntArray(sparseSALength);
            for (int i = 0; i < sparseSA.size(); ++i) {
                sparseSA.set(i, 1L | d.readInt());
            }
            return new SparseSuffixArray(sparseSA, N, L);
        }
        finally {
            if (in != null)
                in.close();
        }
    }

    public static SparseSuffixArray buildFromSuffixArray(LIntArray SA, int L) {
        long sparseSA_length = (SA.size() + L - 1) / L;
        LIntArray sparseSA = new LIntArray(sparseSA_length);
        for (long i = 0; i < sparseSA_length; ++i) {
            sparseSA.set(i, SA.get(i * L));
        }
        return new SparseSuffixArray(sparseSA, SA.size(), L);
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
            return sparseSA.get(pos);

        long cursor = index;
        final long N = fmIndex.textSize();
        for (long j = 1; j <= N; j++) {
            cursor = fmIndex.suffixLink(cursor);
            if (cursor == 0) {
                return j - 1;
            }
            if (cursor % L == 0)
                return sparseSA.get(cursor / L) + j;
        }
        throw new IllegalStateException(String.format("cannot reach here: get(index:%d)", index));

    }
}
