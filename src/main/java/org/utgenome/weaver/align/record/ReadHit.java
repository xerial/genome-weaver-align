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
// ReadHit.java
// Since: 2011/08/25
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.record;

import org.utgenome.weaver.align.Strand;
import org.xerial.lens.SilkLens;

/**
 * Candidate position of read hit
 * 
 * @author leo
 * 
 */
public class ReadHit
{
    public final long    pos;
    public final int     matchLength;
    public final int     diff;
    public final Strand  strand;
    public final boolean isUnique;
    public final ReadHit nextSplit;

    public ReadHit(long pos, int matchLength, int diff, Strand strand, boolean isUnique, ReadHit nextSplit) {
        this.pos = pos;
        this.matchLength = matchLength;
        this.diff = diff;
        this.strand = strand;
        this.isUnique = isUnique;
        this.nextSplit = nextSplit;
    }

    public int getK() {
        if (nextSplit == null)
            return diff;
        else
            return diff + nextSplit.getK();
    }

    public ReadHit addSplit(ReadHit split) {
        return new ReadHit(pos, matchLength, diff, strand, isUnique, split);
    }

    @Override
    public String toString() {
        return SilkLens.toSilk(this);
    }
}
