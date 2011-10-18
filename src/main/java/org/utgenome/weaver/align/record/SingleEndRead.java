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
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.SOLiDColorSequence;

/**
 * 
 * 
 * @author leo
 * 
 */
public class SingleEndRead implements Read
{
    public final String       name;

    /**
     * Read sequence
     */
    public final ACGTSequence seq;

    /**
     * Phread quality value of the read sequence
     */
    public final String       qual;

    public SingleEndRead(String name, String seq, String qual) {
        this.name = name;
        this.seq = new ACGTSequence(seq);
        this.qual = qual;
    }

    public SingleEndRead(String name, ACGTSequence seq, String qual) {
        this.name = name;
        this.seq = seq;
        this.qual = qual;
    }

    @Override
    public String name() {
        return name;
    }

    @Override
    public String toString() {
        return String.format("%s\t%s\t%s", name, seq, qual);
    }

    public static SingleEndRead createFrom(FASTASequence seq) {
        return new SingleEndRead(seq.getSequenceName(), seq.getSequence(), null);
    }

    public static SingleEndRead createFrom(FastqRead seq) {
        return new SingleEndRead(seq.seqname, seq.seq, seq.qual);
    }

    @Override
    public int getNumReadFragment() {
        return 1;
    }

    @Override
    public ACGTSequence getRead(int index) {
        if (index != 0)
            throw new ArrayIndexOutOfBoundsException(index);
        return seq;
    }

    @Override
    public String getQual(int index) {
        if (index != 0)
            throw new ArrayIndexOutOfBoundsException(index);
        return qual;
    }

    @Override
    public boolean isLetterSpace() {
        return true;
    }

    @Override
    public boolean isColorSpace() {
        return false;
    }

    @Override
    public SOLiDColorSequence getColorRead(int index) {
        return null;
    }
}
