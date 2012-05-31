/*
 * sais.java for sais-java
 * Copyright (c) 2008-2010 Yuta Mori All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.utgenome.weaver.align.sais;

import org.utgenome.weaver.align.LIntArray;
import org.utgenome.weaver.align.LSeq;
import org.xerial.util.log.Logger;

/**
 * Suffix-array construction algorithm based on induced-sorting. This
 * implementation uses SAIS by Yuta Mori
 * 
 * @author leo
 * 
 */
public class LSAIS
{
    public static class IntArray implements LSeq
    {
        private int[] m_A   = null;
        private int   m_pos = 0;

        public IntArray(int[] A, int pos) {
            m_A = A;
            m_pos = pos;
        }

        public long lookup(long i) {
            return m_A[(int) (m_pos + i)];
        }

        public void set(long i, long val) {
            m_A[(int) (m_pos + i)] = (int) val;
        }

        public long increment(long i, long val) {
            return m_A[(int) (m_pos + i)] += (int) val;
        }

        @Override
        public long textSize() {
            return m_A.length - m_pos;
        }

        @Override
        public String toString() {
            StringBuilder w = new StringBuilder();
            w.append("[");
            for (long i = 0; i < textSize(); ++i) {
                if (i != 0)
                    w.append(", ");
                w.append(lookup(i));
            }
            w.append("]");
            return w.toString();
        }

    }

    private static class SubIntArray implements LSeq
    {
        private LSeq m_A   = null;
        private long m_pos = 0;

        SubIntArray(LSeq A, long pos) {
            m_A = A;
            m_pos = pos;
        }

        public long lookup(long i) {
            return m_A.lookup(m_pos + i);
        }

        public void set(long i, long val) {
            m_A.set(m_pos + i, val);
        }

        public long increment(long i, long val) {
            return m_A.increment(m_pos + i, val);
        }

        @Override
        public long textSize() {
            return m_A.textSize() - m_pos;
        }
    }

    private static class StringArray implements LSeq
    {
        private String m_A   = null;
        private int    m_pos = 0;

        StringArray(String A, int pos) {
            m_A = A;
            m_pos = pos;
        }

        @Override
        public long textSize() {
            return m_A.length() - m_pos;
        }

        @Override
        public long lookup(long i) {
            return (m_A.charAt((int) (m_pos + i)) & 0xffff);
        }

        @Override
        public void set(long i, long val) {

        }

        @Override
        public long increment(long i, long val) {
            return 0;
        }
    }

    /* find the start or end of each bucket */
    private static void getCounts(LSeq T, LSeq C, long n, int k) {
        int i;
        for (i = 0; i < k; ++i) {
            C.set(i, 0);
        }
        for (i = 0; i < n; ++i) {
            C.increment(T.lookup(i), 1);
        }
    }

    private static void getBuckets(LSeq C, LSeq B, int k, boolean end) {
        int i, sum = 0;
        if (end != false) {
            for (i = 0; i < k; ++i) {
                sum += C.lookup(i);
                B.set(i, sum);
            }
        }
        else {
            for (i = 0; i < k; ++i) {
                sum += C.lookup(i);
                B.set(i, sum - C.lookup(i));
            }
        }
    }

    /* sort all type LMS suffixes */
    private static void LMSsort(LSeq T, LSeq SA, LSeq C, LSeq B, long n, int k) {
        long b, i, j;
        long c0, c1;
        /* compute SAl */
        if (C == B) {
            getCounts(T, C, n, k);
        }
        getBuckets(C, B, k, false); /* find starts of buckets */
        j = n - 1;
        b = B.lookup(c1 = T.lookup(j));
        --j;
        SA.set(b++, (T.lookup(j) < c1) ? ~j : j);
        for (i = 0; i < n; ++i) {
            if (0 < (j = SA.lookup(i))) {
                if ((c0 = T.lookup(j)) != c1) {
                    B.set(c1, b);
                    b = B.lookup(c1 = c0);
                }
                --j;
                SA.set(b++, (T.lookup(j) < c1) ? ~j : j);
                SA.set(i, 0);
            }
            else if (j < 0) {
                SA.set(i, ~j);
            }
        }
        /* compute SAs */
        if (C == B) {
            getCounts(T, C, n, k);
        }
        getBuckets(C, B, k, true); /* find ends of buckets */
        for (i = n - 1, b = B.lookup(c1 = 0); 0 <= i; --i) {
            if (0 < (j = SA.lookup(i))) {
                if ((c0 = T.lookup(j)) != c1) {
                    B.set(c1, b);
                    b = B.lookup(c1 = c0);
                }
                --j;
                SA.set(--b, (T.lookup(j) > c1) ? ~(j + 1) : j);
                SA.set(i, 0);
            }
        }
    }

    private static int LMSpostproc(LSeq T, LSeq SA, long n, long m) {
        long i, j, p, q, plen, qlen;
        int name;
        long c0, c1;
        boolean diff;

        /* compact all the sorted substrings into the first m items of SA
            2*m must be not larger than n (proveable) */
        for (i = 0; (p = SA.lookup(i)) < 0; ++i) {
            SA.set(i, ~p);
        }
        if (i < m) {
            for (j = i, ++i;; ++i) {
                if ((p = SA.lookup(i)) < 0) {
                    SA.set(j++, ~p);
                    SA.set(i, 0);
                    if (j == m) {
                        break;
                    }
                }
            }
        }

        /* store the length of all substrings */
        i = n - 1;
        j = n - 1;
        c0 = T.lookup(n - 1);
        do {
            c1 = c0;
        }
        while ((0 <= --i) && ((c0 = T.lookup(i)) >= c1));
        for (; 0 <= i;) {
            do {
                c1 = c0;
            }
            while ((0 <= --i) && ((c0 = T.lookup(i)) <= c1));
            if (0 <= i) {
                SA.set(m + ((i + 1) >> 1), j - i);
                j = i + 1;
                do {
                    c1 = c0;
                }
                while ((0 <= --i) && ((c0 = T.lookup(i)) >= c1));
            }
        }

        /* find the lexicographic names of all substrings */
        for (i = 0, name = 0, q = n, qlen = 0; i < m; ++i) {
            p = SA.lookup(i);
            plen = SA.lookup(m + (p >> 1));
            diff = true;
            if ((plen == qlen) && ((q + plen) < n)) {
                for (j = 0; (j < plen) && (T.lookup(p + j) == T.lookup(q + j)); ++j) {}
                if (j == plen) {
                    diff = false;
                }
            }
            if (diff != false) {
                ++name;
                q = p;
                qlen = plen;
            }
            SA.set(m + (p >> 1), name);
        }

        return name;
    }

    /* compute SA and BWT */
    private static void induceSA(LSeq T, LSeq SA, LSeq C, LSeq B, long n, int k) {
        long b, i, j;
        long c0, c1;
        /* compute SAl */
        if (C == B) {
            getCounts(T, C, n, k);
        }
        getBuckets(C, B, k, false); /* find starts of buckets */
        j = n - 1;
        b = B.lookup(c1 = T.lookup(j));
        SA.set(b++, ((0 < j) && (T.lookup(j - 1) < c1)) ? ~j : j);
        for (i = 0; i < n; ++i) {
            j = SA.lookup(i);
            SA.set(i, ~j);
            if (0 < j) {
                if ((c0 = T.lookup(--j)) != c1) {
                    B.set(c1, b);
                    b = B.lookup(c1 = c0);
                }
                SA.set(b++, ((0 < j) && (T.lookup(j - 1) < c1)) ? ~j : j);
            }
        }
        /* compute SAs */
        if (C == B) {
            getCounts(T, C, n, k);
        }
        getBuckets(C, B, k, true); /* find ends of buckets */
        for (i = n - 1, b = B.lookup(c1 = 0); 0 <= i; --i) {
            if (0 < (j = SA.lookup(i))) {
                if ((c0 = T.lookup(--j)) != c1) {
                    B.set(c1, b);
                    b = B.lookup(c1 = c0);
                }
                SA.set(--b, ((j == 0) || (T.lookup(j - 1) > c1)) ? ~j : j);
            }
            else {
                SA.set(i, ~j);
            }
        }
    }

    private static Logger _logger = Logger.getLogger(LSAIS.class);

    /* find the suffix array SA of T[0..n-1] in {0..k-1}^n
       use a working space (excluding T and SA) of at most 2n+O(1) for a constant alphabet */
    private static long SA_IS(LSeq T, LSeq SA, long fs, long n, int k) {
        LSeq C, B, RA;
        long b, m, i, j, c, p, q, pidx = 0, newfs;
        int name;
        long c0, c1;
        int flags = 0;

        if (k <= 256) {
            C = new IntArray(new int[k], 0);
            if (k <= fs) {
                B = new SubIntArray(SA, n + fs - k);
                flags = 1;
            }
            else {
                B = new IntArray(new int[k], 0);
                flags = 3;
            }
        }
        else if (k <= fs) {
            C = new SubIntArray(SA, n + fs - k);
            if (k <= (fs - k)) {
                B = new SubIntArray(SA, n + fs - k * 2);
                flags = 0;
            }
            else if (k <= 1024) {
                B = new IntArray(new int[k], 0);
                flags = 2;
            }
            else {
                B = C;
                flags = 8;
            }
        }
        else {
            C = B = new IntArray(new int[k], 0);
            flags = 4 | 8;
        }

        /* stage 1: reduce the problem by at least 1/2
           sort all the LMS-substrings */
        getCounts(T, C, n, k);
        getBuckets(C, B, k, true); /* find ends of buckets */
        for (i = 0; i < n; ++i) {
            SA.set(i, 0);
        }
        b = -1;
        i = n - 1;
        j = n;
        m = 0;
        c0 = T.lookup(n - 1);
        do {
            c1 = c0;
        }
        while ((0 <= --i) && ((c0 = T.lookup(i)) >= c1));
        for (; 0 <= i;) {
            do {
                c1 = c0;
            }
            while ((0 <= --i) && ((c0 = T.lookup(i)) <= c1));
            if (0 <= i) {
                if (0 <= b) {
                    SA.set(b, j);
                }
                b = B.increment(c1, -1);
                j = i;
                ++m;
                do {
                    c1 = c0;
                }
                while ((0 <= --i) && ((c0 = T.lookup(i)) >= c1));
            }
        }
        if (1 < m) {
            LMSsort(T, SA, C, B, n, k);
            name = LMSpostproc(T, SA, n, m);
        }
        else if (m == 1) {
            SA.set(b, j + 1);
            name = 1;
        }
        else {
            name = 0;
        }

        /* stage 2: solve the reduced problem
           recurse if names are not yet unique */
        if (name < m) {
            if ((flags & 4) != 0) {
                C = null;
                B = null;
            }
            if ((flags & 2) != 0) {
                B = null;
            }
            newfs = (n + fs) - (m * 2);
            if ((flags & (1 | 4 | 8)) == 0) {
                if ((k + name) <= newfs) {
                    newfs -= k;
                }
                else {
                    flags |= 8;
                }
            }
            for (i = m + (n >> 1) - 1, j = m * 2 + newfs - 1; m <= i; --i) {
                if (SA.lookup(i) != 0) {
                    SA.set(j--, SA.lookup(i) - 1);
                }
            }
            RA = new SubIntArray(SA, m + newfs);
            SA_IS(RA, SA, newfs, m, name);
            RA = null;

            i = n - 1;
            j = m * 2 - 1;
            c0 = T.lookup(n - 1);
            do {
                c1 = c0;
            }
            while ((0 <= --i) && ((c0 = T.lookup(i)) >= c1));
            for (; 0 <= i;) {
                do {
                    c1 = c0;
                }
                while ((0 <= --i) && ((c0 = T.lookup(i)) <= c1));
                if (0 <= i) {
                    SA.set(j--, i + 1);
                    do {
                        c1 = c0;
                    }
                    while ((0 <= --i) && ((c0 = T.lookup(i)) >= c1));
                }
            }

            for (i = 0; i < m; ++i) {
                SA.set(i, SA.lookup(m + SA.lookup(i)));
            }
            if ((flags & 4) != 0) {
                C = B = new IntArray(new int[k], 0);
            }
            if ((flags & 2) != 0) {
                B = new IntArray(new int[k], 0);
            }
        }

        /* stage 3: induce the result for the original problem */
        if ((flags & 8) != 0) {
            getCounts(T, C, n, k);
        }
        /* put all left-most S characters into their buckets */
        if (1 < m) {
            getBuckets(C, B, k, true); /* find ends of buckets */
            i = m - 1;
            j = n;
            p = SA.lookup(m - 1);
            c1 = T.lookup(p);
            do {
                q = B.lookup(c0 = c1);
                while (q < j) {
                    SA.set(--j, 0);
                }
                do {
                    SA.set(--j, p);
                    if (--i < 0) {
                        break;
                    }
                    p = SA.lookup(i);
                }
                while ((c1 = T.lookup(p)) == c0);
            }
            while (0 <= i);
            while (0 < j) {
                SA.set(--j, 0);
            }
        }

        induceSA(T, SA, C, B, n, k);

        C = null;
        B = null;
        return pidx;
    }

    /** Suffix Sorting **/
    public static long suffixsort(LSeq T, LSeq SA, int k) {
        if (T == null || SA == null)
            throw new NullPointerException();
        final int n = (int) T.textSize();
        if (SA.textSize() < n)
            throw new IllegalArgumentException("The suffix array container (SA) size is smaller than the input");

        if (n <= 1) {
            if (n == 1) {
                SA.set(0, 0);
            }
            return 0;
        }
        return SA_IS(T, SA, 0, n, k);
    }

    public static long suffixsort(LSeq T, int[] SA, int k) {
        if (T == null || SA == null)
            throw new NullPointerException();
        final int n = (int) T.textSize();
        if (SA.length < n)
            throw new IllegalArgumentException("The suffix array container (SA) size is smaller than the input");

        IntArray iSA = new IntArray(SA, 0);

        if (n <= 1) {
            if (n == 1) {
                iSA.set(0, 0);
            }
            return 0;
        }
        return SA_IS(T, iSA, 0, n, k);
    }

    /* String */
    public static long suffixsort(String T, LIntArray SA) {
        if (T == null || SA == null)
            throw new NullPointerException();
        final int n = T.length();
        if (SA.textSize() < n)
            throw new IllegalArgumentException("The suffix array container (SA) size is smaller than the input");

        if (n <= 1) {
            if (n == 1) {
                SA.set(0, 0);
            }
            return 0;
        }
        return SA_IS(new StringArray(T, 0), SA, 0, n, 65536);
    }

}
