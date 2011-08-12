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

    private int[][]            occTable;
    private final ACGTSequence seq;
    private final int          W;
    private final int          K;

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

        _logger.trace("preparing occurrence count table...");
        this.K = ACGT.values().length;
        final int numRows = (int) ((seq.textSize() + W) / W);

        occTable = new int[numRows][K];
        for (int k = 0; k < K; ++k) {
            occTable[0][k] = 0;
        }
        for (int i = 1; i < numRows; ++i) {
            long start = (i - 1) * W;
            long end = Math.min(start + W, seq.textSize());
            for (ACGT each : ACGT.values()) {
                occTable[i][each.code] = occTable[i - 1][each.code] + (int) seq.fastCount(each, start, end);
            }
        }
        _logger.trace("done.");
    }

    /**
     * Get the character occurrence counts of ACGTN in seq[0..index)
     * 
     * @param index
     * @return
     */
    public int[] getOccACGT(long index) {
        if (index > seq.textSize())
            index = seq.textSize();
        int blockPos = (int) (index / W);

        int[] occ = new int[K];
        for (int i = 0; i < K; ++i) {
            occ[i] = occTable[blockPos][i] + (int) seq.fastCount(ACGT.decode((byte) i), blockPos * W, index);
        }
        return occ;
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
        if (index > seq.textSize())
            index = seq.textSize();
        int blockPos = (int) (index / W);
        // Look up the occurrence count table. 
        // And also count the characters using the original sequence
        int occ = occTable[blockPos][ch.code] + (int) seq.fastCount(ch, blockPos * W, index);
        return occ;
    }

}
