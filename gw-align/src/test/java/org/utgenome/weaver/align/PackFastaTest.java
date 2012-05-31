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
// PackFastaTest.java
// Since: 2011/07/13
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.io.File;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.utgenome.util.TestHelper;
import org.utgenome.weaver.GenomeWeaver;
import org.xerial.util.FileUtil;
import org.xerial.util.log.Logger;

public class PackFastaTest
{
    private static Logger _logger = Logger.getLogger(PackFastaTest.class);

    File                  tmpDir;

    @Before
    public void setUp() throws Exception {
        tmpDir = TestHelper.createTempDir();
    }

    @After
    public void tearDown() throws Exception {
        if (tmpDir != null && tmpDir.exists()) {
            FileUtil.rmdir(tmpDir);
        }
    }

    @Test
    public void pack() throws Exception {

        File fastaArchive = TestHelper.createTempFileFrom(PackFastaTest.class, "test.fa", new File(tmpDir, "test.fa"));
        GenomeWeaver.execute(String.format("PackFasta %s", fastaArchive));

        BWTFiles fdb = new BWTFiles(fastaArchive.getPath(), Strand.FORWARD);
        BWTFiles rdb = new BWTFiles(fastaArchive.getPath(), Strand.REVERSE);
        ACGTSequence seq = ACGTSequence.loadFrom(fdb.pac());
        _logger.info(seq);
        ACGTSequence rseq = ACGTSequence.loadFrom(rdb.pac());
        _logger.info(rseq);
    }

    @Test
    public void packTarGZ() throws Exception {

        File fastaArchive = TestHelper.createTempFileFrom(BWTransformTest.class, "sample-archive.fa.tar.gz", new File(
                tmpDir, "sample.fa.tar.gz"));
        GenomeWeaver.execute(String.format("PackFasta %s", fastaArchive));

        BWTFiles fdb = new BWTFiles(fastaArchive.getPath(), Strand.FORWARD);
        BWTFiles rdb = new BWTFiles(fastaArchive.getPath(), Strand.REVERSE);
        ACGTSequence seq = ACGTSequence.loadFrom(fdb.pac());
        _logger.info(seq);
        ACGTSequence rseq = ACGTSequence.loadFrom(rdb.pac());
        _logger.info(rseq);
    }

}
