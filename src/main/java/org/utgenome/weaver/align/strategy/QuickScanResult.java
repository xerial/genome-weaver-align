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
// QuickScanResult.java
// Since: 2011/08/13
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.strategy;

import org.utgenome.weaver.align.ACGT;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.BitVector;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.Range;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.SuffixInterval;
import org.xerial.lens.SilkLens;

/**
 * For quick lookup of the mismatch locations in a read sequence
 * 
 * @author leo
 * 
 */
public class QuickScanResult
{
    public final SuffixInterval si;
    public final BitVector      breakPoint;
    public final int            numMismatches;
    public final Range          longestMatch;
    public final SuffixInterval longestMatchSi;

    public QuickScanResult(SuffixInterval si, BitVector breakPoint, int numMismatches, Range longestMatch,
            SuffixInterval longestMatchSi) {
        this.si = si;
        this.breakPoint = breakPoint;
        this.numMismatches = numMismatches;
        this.longestMatch = longestMatch;
        this.longestMatchSi = longestMatchSi;
    }

    @Override
    public String toString() {
        return SilkLens.toSilk(this);
    }

    public static QuickScanResult scanMismatchLocations(FMIndexOnGenome fmIndex, ACGTSequence query, Strand strand) {
        int qLen = (int) query.textSize();
        int numMismatches = 0;
        BitVector breakPoint = new BitVector(qLen);
        SuffixInterval si = fmIndex.wholeSARange();
        int longestMatchLength = 0;
        int mark = 0;
        Range longestMatch = null;
        SuffixInterval longestMatchSi = null;
        int i = 0;
        for (; i < qLen; ++i) {
            ACGT ch = query.getACGT(i);
            si = fmIndex.forwardSearch(strand, ch, si);
            if (si.isEmpty()) {
                breakPoint.set(i, true);
                numMismatches++;
                if (longestMatch == null || longestMatch.length() < (i - mark)) {
                    longestMatch = new Range(mark, i);
                    longestMatchSi = si;
                }
                si = fmIndex.wholeSARange();
                mark = i + 1;
            }
        }
        if (longestMatch == null || longestMatch.length() < (i - mark)) {
            longestMatch = new Range(mark, i);
        }

        return new QuickScanResult(si, breakPoint, numMismatches, longestMatch, longestMatchSi);
    }
}
