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
// BurrowsWheelerAlignmentTest.java
// Since: 2011/02/10
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
import org.utgenome.weaver.align.record.AlignmentRecord;
import org.xerial.lens.SilkLens;
import org.xerial.util.ObjectHandlerBase;
import org.xerial.util.log.Logger;

public class BWAlignTest
{
    private static Logger _logger = Logger.getLogger(BWAlignTest.class);

    File                  tmpDir;

    @Before
    public void setUp() throws Exception {
        tmpDir = new File("target", "bwa");
        if (!tmpDir.exists())
            tmpDir.mkdirs();
    }

    @After
    public void tearDown() throws Exception {
        if (tmpDir != null && tmpDir.exists()) {
            //FileUtil.rmdir(tmpDir);
        }
    }

    @Test
    public void align() throws Exception {

        File fastaArchive = TestHelper.createTempFileFrom(BWTTest.class, "test.fa", new File(tmpDir, "test.fa"));
        GenomeWeaver.execute(String.format("bwt %s", fastaArchive));
        GenomeWeaver.execute(String.format("align %s -q TATAA", fastaArchive.getPath()));

    }

    @Test
    public void align2() throws Exception {
        File fastaArchive = TestHelper.createTempFileFrom(BWTTest.class, "test2.fa", new File(tmpDir, "test.fa"));
        GenomeWeaver.execute(String.format("bwt %s", fastaArchive));
        GenomeWeaver.execute(String.format("align %s -q ATCTCATGGGA", fastaArchive.getPath()));
    }

    @Test
    public void align3() throws Exception {
        File fastaArchive = TestHelper.createTempFileFrom(BWTTest.class, "test2.fa", new File(tmpDir, "test2.fa"));
        GenomeWeaver.execute(String.format("bwt %s", fastaArchive));

        BWAlign.query(fastaArchive.getPath(), "TAAAGTAT", new ObjectHandlerBase<AlignmentRecord>() {
            @Override
            public void handle(AlignmentRecord input) throws Exception {
                _logger.info(SilkLens.toSilk(input));
                assertEquals("seq2", input.chr);
                assertEquals(Strand.REVERSE, input.strand);
                assertEquals(9, input.start);
                assertEquals(17, input.end); // 1-origin
            }
        });

        BWAlign.query(fastaArchive.getPath(), "ATACTTTA", new ObjectHandlerBase<AlignmentRecord>() {
            @Override
            public void handle(AlignmentRecord input) throws Exception {
                _logger.info(SilkLens.toSilk(input));
                assertEquals("seq2", input.chr);
                assertEquals(Strand.FORWARD, input.strand);
                assertEquals(9, input.start); // 1-origin
                assertEquals(17, input.end); // 1-origin
            }
        });

    }

    @Test
    public void sample() throws Exception {
        File fastaArchive = TestHelper.createTempFileFrom(BWTTest.class, "sample.fa", new File(tmpDir, "sample.fa"));
        GenomeWeaver.execute(String.format("bwt %s", fastaArchive));

        FASTAPullParser fa = new FASTAPullParser(fastaArchive);
        final FASTASequence seq = fa.nextSequence();
        fa.close();

        BWAlign.query(fastaArchive.getPath(), "TTTCAG", new ObjectHandlerBase<AlignmentRecord>() {
            @Override
            public void handle(AlignmentRecord input) throws Exception {
                _logger.debug(SilkLens.toSilk(input));
                String s = seq.getSequence().substring(input.start - 1, input.end - 1);
                assertEquals(String.format("strand:%s query:%s ref:%s", input.strand, input.querySeq, s), s,
                        input.querySeq);
            }
        });

    }

    @Test
    public void sample2() throws Exception {
        File fastaArchive = TestHelper.createTempFileFrom(BWTTest.class, "sample2.fa", new File(tmpDir, "sample2.fa"));
        GenomeWeaver.execute(String.format("bwt %s", fastaArchive));

        FASTAPullParser fa = new FASTAPullParser(fastaArchive);
        final FASTASequence seq = fa.nextSequence();
        fa.close();

        BWAlign.query(fastaArchive.getPath(), "TTTCAG", new ObjectHandlerBase<AlignmentRecord>() {
            @Override
            public void handle(AlignmentRecord input) throws Exception {
                _logger.info(SilkLens.toSilk(input));
                String s = seq.getSequence().substring(input.start - 1, input.end - 1);
                assertEquals(String.format("strand:%s query:%s ref:%s", input.strand, input.querySeq, s), s,
                        input.querySeq);
            }
        });

    }

}
