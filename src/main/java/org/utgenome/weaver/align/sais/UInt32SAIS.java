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
// Uint32SAIS.java
// Since: 2011/02/17
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.sais;

import org.utgenome.weaver.align.LSeq;
import org.utgenome.weaver.align.RSBitVector;

/**
 * Suffix-array (SA) construction algorithm based on Induced Sorting (IS)
 * 
 * @author leo
 * 
 */
public class UInt32SAIS
{
    private final LSeq           T;
    private final long           N;
    private final int            K;
    private final long[]         bucket;
    private RSBitVector          typeLS;

    private final static boolean LType = false;
    private final static boolean SType = true;

    static class ArrayWrap implements LSeq
    {

        private final LSeq seq;
        private final long offset;
        private final long length;

        public ArrayWrap(LSeq original, long offset) {
            this.seq = original;
            this.offset = offset;
            this.length = original.textSize() - offset;
        }

        public ArrayWrap(LSeq original, long offset, long length) {
            this.seq = original;
            this.offset = offset;
            this.length = length;
        }

        @Override
        public long lookup(long index) {
            return seq.lookup(index + offset);
        }

        @Override
        public long textSize() {
            return length;
        }

        @Override
        public void set(long index, long value) {
            seq.set(index + offset, value);
        }

        @Override
        public long update(long index, long value) {
            return seq.update(index + offset, value);
        }

        @Override
        public String toString() {
            StringBuilder w = new StringBuilder();
            w.append("[");
            for (long i = 0; i < length; ++i) {
                if (i != 0)
                    w.append(", ");
                w.append(lookup(i));
            }
            w.append("]");
            return w.toString();
        }
    }

    public UInt32SAIS(LSeq T, int K) {
        this.T = T;
        this.N = T.textSize();
        this.K = K;
        this.bucket = new long[K + 1];
        typeLS = new RSBitVector(N);
    }

    public static UInt32Array SAIS(LSeq T, int K) {
        UInt32Array SA = new UInt32Array(T.textSize());
        SAIS(T, SA, K);
        return SA;
    }

    public static LSeq SAIS(LSeq T, LSeq SA, int K) {
        UInt32SAIS sais = new UInt32SAIS(T, K);
        sais.SAIS(SA);
        return SA;
    }

    public void SAIS(LSeq SA) {

        // initialize the suffix array
        for (long i = 0; i < N; ++i)
            SA.set(i, 0);

        typeLS.setBit(SType, N - 1); // the sentinel 
        typeLS.setBit(LType, N - 2);

        // set the type of each character
        for (long i = N - 2; i > 0; --i) {
            long x = T.lookup(i);
            long y = T.lookup(i - 1);
            if (x < y)
                typeLS.setBit(LType, i - 1);
            else if (x > y)
                typeLS.setBit(SType, i - 1);
            else
                typeLS.setBit(typeLS.get(i), i - 1);
        }

        // Step 1: reduce the problem by at least 1/2 
        // sort all the S-substrings

        findEndOfBuckets();

        // Find LMS characters
        for (int i = 1; i < N; ++i) {
            if (isLMS(i))
                SA.set(--bucket[(int) T.lookup(i)], i);
        }

        induceSA_left(SA);
        induceSA_right(SA);

        // Compact all the sorted substrings into the first M items of SA
        // 2*M must be not larger than N 
        int numLMS = 0;
        for (long i = 0; i < N; ++i) {
            if (isLMS(SA.lookup(i)))
                SA.set(numLMS++, SA.lookup(i));
        }

        // Initialize the name array buffer
        for (long i = numLMS; i < N; ++i)
            SA.set(i, 0);

        // Find the lexicographic names of the LMS substrings
        int name = 0;
        long prev = -1;
        long qlen = 0;
        for (long i = 0; i < numLMS; ++i) {
            final long pos = SA.lookup(i);
            //final long plen = SA.lookup(numLMS + (pos >> 1));
            boolean diff = false;

            for (long d = 0; d < N; ++d) {
                if (prev == -1 || T.lookup(pos + d) != T.lookup(prev + d)
                        || typeLS.get(pos + d) != typeLS.get(prev + d)) {
                    diff = true;
                    break;
                }
                else if (d > 0 && (isLMS(pos + d) || isLMS(prev + d)))
                    break;
            }

            if (diff) {
                ++name;
                prev = pos;
            }

            SA.set(numLMS + (pos >> 1), name - 1);
        }

        for (long i = N - 1, j = N - 1; i >= numLMS; --i) {
            if (SA.lookup(i) != 0)
                SA.set(j--, SA.lookup(i));
        }

        // Step 2: solve the reduced problem
        // Build SA1 
        LSeq SA1 = new ArrayWrap(SA, 0, N - numLMS);
        LSeq T1 = new ArrayWrap(SA, N - numLMS, numLMS);
        if (name < numLMS) {
            new UInt32SAIS(T1, name - 1).SAIS(SA1);
        }
        else {
            // Generate the suffix array of inputS1 directory.
            for (long i = 0; i < numLMS; i++)
                SA1.set(T1.lookup(i), i);
        }

        // Step 3: Induce the result for the original problem
        findEndOfBuckets();
        for (long i = 1, j = 0; i < N; ++i) {
            if (isLMS(i))
                T1.set(j++, i); // get p1
        }
        for (long i = 0; i < numLMS; ++i) {
            SA1.set(i, T1.lookup(SA1.lookup(i)));
        }
        // init SA[N1 .. N-1]
        for (long i = numLMS; i < N; ++i) {
            SA.set(i, 0);
        }
        for (long i = numLMS - 1; i >= 0; --i) {
            long j = SA.lookup(i);
            SA.set(i, 0);
            SA.set(--bucket[(int) T.lookup(j)], j);
        }
        induceSA_left(SA);
        induceSA_right(SA);

    }

    private void findStartOfBuckets() {
        initBuckets();
        // compute the start of the buckets
        int sum = 0;
        for (int i = 0; i <= K; ++i) {
            sum += bucket[i];
            bucket[i] = sum - bucket[i];
        }
    }

    private void findEndOfBuckets() {
        initBuckets();
        // compute the end of the buckets
        long sum = 0;
        for (int i = 0; i <= K; ++i) {
            sum += bucket[i];
            bucket[i] = sum;
        }
    }

    private void initBuckets() {
        // initialize buckets
        for (int i = 0; i <= K; ++i) {
            bucket[i] = 0;
        }
        // compute the size of each bucket
        for (int i = 0; i < N; ++i) {
            ++bucket[(int) T.lookup(i)];
        }
    }

    boolean isLMS(long pos) {
        return typeLS.get(pos) && !typeLS.get(pos - 1);
    }

    private void induceSA_left(LSeq SA) {
        findStartOfBuckets();
        long j;
        for (long i = 0; i < N; ++i) {
            j = SA.lookup(i) - 1;
            if (j >= 0 && !typeLS.get(j))
                SA.set(bucket[(int) T.lookup(j)]++, j);
        }
    }

    private void induceSA_right(LSeq SA) {
        findEndOfBuckets();
        long j;
        for (long i = N - 1; i >= 0; --i) {
            j = SA.lookup(i) - 1;
            if (j >= 0 && typeLS.get(j))
                SA.set(--bucket[(int) T.lookup(j)], j);
        }
    }

}
