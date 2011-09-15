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
    public final BitVector            chunkWithMismatch;

    public PrefixScan(Strand strand, List<SuffixInterval> si, BitVector matchedChunkFlag) {
        this.strand = strand;
        this.si = si;
        this.chunkWithMismatch = matchedChunkFlag;
    }

    public List<Long> hitCount() {
        List<Long> count = new ArrayList<Long>();
        for (SuffixInterval each : si) {
            if (each == null)
                count.add(0L);
            else
                count.add(each.range());
        }
        return count;
    }

    @Override
    public String toString() {
        return String.format("%s [%s]", chunkWithMismatch.toStringReverse(), StringUtil.join(hitCount(), ", "));
    }

    public static PrefixScan scanRead(FMIndexOnGenome fmIndex, ACGTSequence query, Strand strand, StaircaseFilter filter) {

        final int s = filter.getNumChunks();
        BitVector chunkWithMismatch = new BitVector(s);
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
                    chunkWithMismatch.set(c);
                    continue chunk_loop;
                }
            }
            siOfChunks.add(si);

        }

        return new PrefixScan(strand, siOfChunks, chunkWithMismatch);
    }
}
