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

import org.utgenome.weaver.align.ACGT;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.BitVector;
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

    // automaton
    private BitVector[]           patternMask;
    private BitVector[]           automaton;
    private BitVector[]           stairMask;

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

        // prefix length of each chunk of the filter
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

        patternMask = new BitVector[ACGT.values().length];
        for (int i = 0; i < patternMask.length; ++i)
            patternMask[i] = new BitVector(m);

        // Set bit flag where the character appears.
        for (BitVector each : patternMask) {
            each.clear();
        }
        for (int i = 0; i < m; ++i) {
            ACGT ch = query.getACGT(i);
            patternMask[ch.code].set(i);
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

    /**
     * <pre>
     *   *---*---*---*
     *   | \ | \ | \ |
     *   *---*---*---*
     *   | \ | \ | \ |
     *   *---*---*---*
     * 
     * </pre>
     * 
     * 
     * 
     * @param strand
     * @param out
     * @throws Exception
     */
    public void match(ObjectHandler<Candidate> out) throws Exception {

        // row-wise simulation of NFA 
        for (int row = 0; row <= k; ++row) {
            if (row == k) {
                // No more difference is allowed, so try to find a strong match
                SuffixInterval si = exactMatch(chunkStart[row], chunkStart[row + 1]);
                if (!si.isValidRange())
                    continue;
                out.handle(new Candidate(si, k, 2 * k + m));
            }
            else {
                rewind = chunkStart[row] + k;
                simulateNFA0(row, out);
            }
        }

    }

    private void simulateNFA0(int row, ObjectHandler<Candidate> out) throws Exception {

        // Starts with the chunk position  
        final int startChunk = row;
        // First chunk must have an exact match  
        SuffixInterval si = exactMatch(chunkStart[startChunk], chunkStart[startChunk + 1]);
        if (si.isEmpty()) {
            return;
        }

        // Init the automaton
        int height = k - row + 1; // when k=2 and row=0, the automaton height is 3
        automaton = new BitVector[height];
        for (int i = 0; i < automaton.length; ++i)
            automaton[i] = new BitVector(m);

        // Prepare staircase filter
        stairMask = new BitVector[height];
        for (int i = 0; i < stairMask.length; ++i) {
            stairMask[i] = new BitVector(m).not().rshift(chunkStart[startChunk + i]);
        }

        // Init the diagonal of the automaton
        byte q = chunkStart[startChunk + 1];
        automaton[0].set(q);
        for (int i = 1; i < height; ++i) {
            if (stairMask[i].get(q + 1)) {
                q++;
                automaton[i].set(q);
            }
            else
                break;
        }

        dfsSuffix(0, si, out);

    }

    private void dfsSuffix(int step, SuffixInterval si, ObjectHandler<Candidate> out) throws Exception {
        for (ACGT ch : ACGT.exceptN) {
            SuffixInterval nextSi = fmIndex.forwardSearch(strand, ch, si);
            if (nextSi.isEmpty())
                continue;

            ++rewind;
            DFSState state = activate(step, ch);
            switch (state.type) {
            case Finished:
                out.handle(new Candidate(nextSi, state.rc, rewind));
                break;
            case NoMatch:
                break;
            case Match:
                if (nextSi.range() == 1 && (m - step) < 8) {
                    // unique hit
                    out.handle(new Candidate(nextSi, 0, rewind));
                }
                else {
                    dfsSuffix(step + 1, nextSi, out);
                }
                break;
            }
            --rewind;
        }
    }

    private enum State {
        Finished, NoMatch, Match
    }

    private static class DFSState
    {
        public final State type;
        public final int   rc;

        public DFSState(State state, int rc) {
            this.type = state;
            this.rc = rc;
        }
    }

    private DFSState activate(int step, ACGT ch) {

        return null;
    }

}
