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
// PrefixScan.java
// Since: 2011/09/15
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.strategy;

import java.util.ArrayList;
import java.util.List;

import org.utgenome.weaver.align.ACGT;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.BitVector;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.SuffixInterval;
import org.xerial.util.StringUtil;

public class PrefixScan
{
    public final Strand               strand;
    public final List<SuffixInterval> si;
    public final BitVector            matchedChunkFlag;

    public PrefixScan(Strand strand, List<SuffixInterval> si, BitVector matchedChunkFlag) {
        this.strand = strand;
        this.si = si;
        this.matchedChunkFlag = matchedChunkFlag;
    }

    @Override
    public String toString() {
        return String.format("%s[%s]", matchedChunkFlag, StringUtil.join(si, ", "));
    }

    public static PrefixScan scanRead(FMIndexOnGenome fmIndex, ACGTSequence query, Strand strand, StaircaseFilter filter) {

        final int s = filter.getNumChunks();
        BitVector matchedChunk = new BitVector(s);
        ArrayList<SuffixInterval> siOfChunks = new ArrayList<SuffixInterval>(s);

        // for each chunk
        chunk_loop: for (int c = 0; c < s; ++c) {
            SuffixInterval si = fmIndex.wholeSARange();
            int cursor = filter.getChunkStart(c);
            final int chunkEnd = cursor + filter.getChunkSize(c);
            for (int x = cursor; x < chunkEnd; ++x) {
                ACGT ch = query.getACGT(x);

                si = fmIndex.forwardSearch(strand, ch, si);
                if (si.isEmpty()) {
                    siOfChunks.add(null);
                    continue chunk_loop;
                }
            }
            siOfChunks.add(si);
            matchedChunk.set(c);
        }

        return new PrefixScan(strand, siOfChunks, matchedChunk);
    }
}
