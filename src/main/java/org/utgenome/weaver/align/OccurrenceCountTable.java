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
// OccurrenceTable.java
// Since: 2011/02/10
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.util.ArrayList;

import org.xerial.util.log.Logger;

/**
 * Occurrence count table for FM-index
 * 
 * @author leo
 * 
 */
public class OccurrenceCountTable
{
    private static Logger      _logger = Logger.getLogger(OccurrenceCountTable.class);

    private ArrayList<int[]>   occTable;
    private final ACGTSequence seq;
    private final int          W;

    /**
     * Create a character occurrence count table. This class saves the
     * occurrence counts of each character for every W-bp block.
     * 
     * @param seq
     * @param windowSize
     */
    public OccurrenceCountTable(ACGTSequence seq, int windowSize) {
        this.seq = seq;
        this.W = windowSize;

        _logger.debug("preparing occurrence count table...");
        final int K = ACGT.values().length;
        final int tableSize = (int) ((seq.textSize() + W) / W);

        occTable = new ArrayList<int[]>(K);
        for (int k = 0; k < K; ++k) {
            int[] occ = new int[tableSize];
            for (int i = 0; i < occ.length; ++i)
                occ[i] = 0;
            occTable.add(occ);
        }

        for (int i = 0; i < seq.textSize(); ++i) {
            ACGT base = seq.getACGT(i);
            int codeIndex = base.code;
            final int blockIndex = i / W;

            if (i % W == 0 && blockIndex > 0) {
                for (ACGT each : ACGT.values()) {
                    int[] occ = occTable.get(each.code);
                    occ[blockIndex] = occ[blockIndex - 1];
                }
            }
            occTable.get(codeIndex)[blockIndex]++;
        }
        _logger.trace("done.");
    }

    /**
     * Get the occurrence count of the specified character contained in
     * seq[0..index)
     * 
     * @param ch
     * @param index
     * @return
     */
    public int getOcc(ACGT ch, long index) {
        long blockPos = index / W;
        if (blockPos > Integer.MAX_VALUE) {
            // index / W must be smaller than Integer.MIN 
            throw new IllegalStateException("Occ table size cannot exceed 2^31-1");
        }
        // Look up the occurrence count table
        int occ = blockPos <= 0 ? 0 : occTable.get(ch.code)[(int) (blockPos - 1)];
        // Count the characters using the original sequence
        final long upperLimit = Math.min(index, seq.textSize());
        occ += seq.fastCount(ch, blockPos * W, upperLimit);
        return occ;
    }

}
