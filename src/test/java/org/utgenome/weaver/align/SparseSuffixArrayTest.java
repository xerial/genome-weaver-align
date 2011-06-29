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
// SparseSuffixArrayTest.java
// Since: 2011/04/27
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import static org.junit.Assert.*;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;

import org.junit.Test;
import org.utgenome.format.fasta.FASTAPullParser;
import org.utgenome.format.fasta.FASTASequence;
import org.utgenome.weaver.align.sais.UInt32Array;
import org.utgenome.weaver.align.sais.UInt32SAIS;
import org.xerial.util.FileResource;

public class SparseSuffixArrayTest
{

    @Test
    public void shortSeq() throws Exception {
        IUPACSequence seq = new IUPACSequence("AACCTATCTATACCCCGGGGAATTATATACGTCCTTATACTTTATCTCATGGGATA", true);
        UInt32Array SA = UInt32SAIS.SAIS(seq, 16);
        SparseSuffixArray ssa = SparseSuffixArray.buildFromSuffixArray(SA, 32);

        ByteArrayOutputStream buf = new ByteArrayOutputStream();
        ssa.saveTo(new DataOutputStream(buf));
        buf.close();

        SparseSuffixArray ssa2 = SparseSuffixArray.loadFrom(new DataInputStream(new ByteArrayInputStream(buf
                .toByteArray())));

        assertEquals(ssa.sparseSA.textSize(), ssa2.sparseSA.textSize());
        for (int i = 0; i < ssa.sparseSA.textSize(); ++i) {
            assertEquals(ssa.sparseSA.lookup(i), ssa2.sparseSA.lookup(i));
        }
    }

    @Test
    public void longSeq() throws Exception {

        FASTAPullParser fa = new FASTAPullParser(FileResource.open(SparseSuffixArrayTest.class, "sample.fa"));
        FASTASequence fa1 = fa.nextSequence();
        fa.close();
        IUPACSequence seq = new IUPACSequence(fa1.getSequence(), true);
        UInt32Array SA = UInt32SAIS.SAIS(seq, 16);

        final int L = 16;
        SparseSuffixArray ssa = SparseSuffixArray.buildFromSuffixArray(SA, L);

        for (int i = 0; i < ssa.sparseSA.textSize(); ++i) {
            assertEquals(SA.lookup(i * L), ssa.sparseSA.lookup(i));
        }

        ByteArrayOutputStream buf = new ByteArrayOutputStream();
        ssa.saveTo(new DataOutputStream(buf));
        buf.close();

        SparseSuffixArray ssa2 = SparseSuffixArray.loadFrom(new DataInputStream(new ByteArrayInputStream(buf
                .toByteArray())));

        assertEquals(ssa.sparseSA.textSize(), ssa2.sparseSA.textSize());
        for (int i = 0; i < ssa.sparseSA.textSize(); ++i) {
            assertEquals(ssa.sparseSA.lookup(i), ssa2.sparseSA.lookup(i));
        }
    }

}
