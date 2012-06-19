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
// BWTransformTest.java
// Since: 2011/04/26
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import static org.junit.Assert.*;

import java.io.File;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.utgenome.format.fasta.FASTAPullParser;
import org.utgenome.format.fasta.FASTASequence;
import org.utgenome.util.TestHelper;
import org.utgenome.weaver.GenomeWeaver;
import org.xerial.util.FileUtil;

public class BWTransformTest
{
    File tmpDir;

    @Before
    public void setUp() throws Exception {
        tmpDir = new File("target", "bwa");
        if (!tmpDir.exists())
            tmpDir.mkdirs();
    }

    @After
    public void tearDown() throws Exception {
        if (tmpDir != null && tmpDir.exists()) {
            FileUtil.rmdir(tmpDir);
        }
    }

    @Test
    public void bwt1() throws Exception {
        File fastaArchive = TestHelper
                .createTempFileFrom(BWTransformTest.class, "test.fa", new File(tmpDir, "test.fa"));
        GenomeWeaver.execute(String.format("bwt %s", fastaArchive));
    }

    @Test
    public void bwt2() throws Exception {
        File fastaArchive = TestHelper.createTempFileFrom(BWTransformTest.class, "sample-archive.fa.tar.gz", new File(
                tmpDir, "sample.fa.tar.gz"));
        GenomeWeaver.execute(String.format("bwt %s", fastaArchive));
    }

    @Test
    public void transform() throws Exception {
        File fastaArchive = TestHelper.createTempFileFrom(BWTransformTest.class, "sample.fa", new File(tmpDir,
                "sample.fa"));
        GenomeWeaver.execute(String.format("bwt %s", fastaArchive));

        FASTAPullParser fa = new FASTAPullParser(fastaArchive);
        FASTASequence orig = fa.nextSequence();

        // Forward
        {
            BWTFiles db = new BWTFiles("target/bwa/sample.fa", Strand.FORWARD);
            ACGTSequence seq = ACGTSequence.loadFrom(db.pac());

            // forward
            for (int i = 0; i < orig.getSequence().length(); ++i) {
                ACGT c = seq.getACGT(i);

                char origBase = Character.toUpperCase(orig.getSequence().charAt(i));
                char iupacBase = Character.toUpperCase(c.name().charAt(0));
                assertEquals(String.format("index %d. orig:%s iupac:%s", i, origBase, c), origBase, iupacBase);
            }
        }

        // Reverse
        {
            BWTFiles db = new BWTFiles("target/bwa/sample.fa", Strand.REVERSE);
            ACGTSequence seq = ACGTSequence.loadFrom(db.pac());
            final int origLen = orig.getSequence().length();

            int offset = 0;

            for (int i = 0; i < orig.getSequence().length(); ++i) {
                ACGT c = seq.getACGT(i + offset);

                char origBase = Character.toUpperCase(orig.getSequence().charAt(origLen - i - 1));
                char iupacBase = Character.toUpperCase(c.name().charAt(0));
                assertEquals(String.format("index:%d, orig:%s, iupac:%s", i, origBase, c.name()), origBase, iupacBase);
            }
        }

    }
}
