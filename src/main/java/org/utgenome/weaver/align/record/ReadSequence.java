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
// ReadSequence.java
// Since: 2011/04/28
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.record;

import org.utgenome.format.fasta.FASTASequence;
import org.utgenome.format.fastq.FastqRead;

/**
 * 
 * 
 * @author leo
 * 
 */
public class ReadSequence implements RawRead
{
    public final String name;

    /**
     * Read sequence
     */
    public final String seq;
    /**
     * Phread quality value of the read sequence
     */
    public final String qual;

    public ReadSequence(String name, String seq, String qual) {
        this.name = name;
        this.seq = seq;
        this.qual = qual;
    }

    @Override
    public String toString() {
        return String.format("%s\t%s\t%s", name, seq, qual);
    }

    public static ReadSequence createFrom(FASTASequence seq) {
        return new ReadSequence(seq.getSequenceName(), seq.getSequence(), null);
    }

    public static ReadSequence createFrom(FastqRead seq) {
        return new ReadSequence(seq.seqname, seq.seq, seq.qual);
    }

}
