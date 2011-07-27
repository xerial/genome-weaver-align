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
// SuffixFilter.java
// Since: 2011/07/27
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.strategy;

import java.util.Arrays;

import org.utgenome.weaver.align.ACGT;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.SuffixInterval;
import org.xerial.util.ObjectHandler;

public class SuffixFilter
{
    private final int             k;
    private final int             m;
    private final FMIndexOnGenome fmIndex;
    private final ACGTSequence    query;
    private final Strand          strand;

    private byte[]                chunkStart;
    private byte[]                chunkLen;
    private int                   rewind;

    public static class Candidate
    {
        public final SuffixInterval si;
        public final int            k;
        public final int            offset;

        public Candidate(SuffixInterval si, int k, int offset) {
            this.si = si;
            this.k = k;
            this.offset = offset;
        }
    }

    public SuffixFilter(int k, FMIndexOnGenome fmIndex, ACGTSequence query, Strand strand) {
        this.k = k;
        this.fmIndex = fmIndex;
        this.query = query;
        this.strand = strand;

        this.m = (int) query.textSize();
        int lastChunkSize = (m - k >= 6) ? m * 2 / (k + 2) : m - k;

        // prefix length of each chunk
        chunkStart = new byte[k + 2];
        byte rest = (byte) (m - lastChunkSize);
        if (k == 0)
            chunkStart[0] = rest;
        else {
            for (int i = 0; i <= k; ++i) {
                chunkStart[i] = (byte) (rest * i / k);
            }
        }
        chunkStart[k + 1] = (byte) m;

        chunkLen = new byte[k + 1];
        for (int i = 0; i <= k; ++i) {
            chunkLen[i] = (byte) (chunkStart[i + 1] - chunkStart[i]);
        }

    }

    public SuffixInterval exactMatch(int start, int end) {
        SuffixInterval si = fmIndex.wholeSARange();
        for (int i = start; i < end; ++i) {
            ACGT ch = query.getACGT(i);
            si = fmIndex.forwardSearch(strand, ch, si);
            if (si.isEmpty())
                break;
        }
        return si;
    }

    public void match(ObjectHandler<Candidate> out) throws Exception {

        for (int i = 0; i < k + 1; ++i) {
            int allowedDiff = k - i;
            if (allowedDiff == 0) {
                // Find strong match
                SuffixInterval si = exactMatch(chunkStart[i], chunkStart[i + 1]);
                if (!si.isValidRange())
                    continue;
                out.handle(new Candidate(si, k, 2 * k + m));
            }
            else {
                rewind = chunkStart[i] + k;
                dfs0(allowedDiff, out);
            }
        }

    }

    private void dfs0(int stairLevel, ObjectHandler<Candidate> out) {
        int chunkPos = k - stairLevel;
        SuffixInterval si = exactMatch(chunkStart[chunkPos], chunkStart[chunkPos + 1]);
        if (si.isEmpty()) {
            return;
        }

        int[] allowedDiff = new int[m + 1];
        Arrays.fill(allowedDiff, 0);
        int nextChunk = chunkPos + 1;
        int e = 0;
        for (int i = chunkStart[chunkPos]; i <= m; ++i) {
            allowedDiff[i] = e;
            if (i == chunkStart[nextChunk]) {
                ++e;
                ++nextChunk;
            }
        }

        int height = stairLevel + 1;
        int maxPrefixLen = m - chunkStart[chunkPos + 1] + stairLevel;

    }
}
