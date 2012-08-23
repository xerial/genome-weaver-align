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
// ConsensusBuilderTest.java
// Since: 2011/10/05
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.assembly;

import static org.junit.Assert.*;

import java.io.File;

import org.junit.Test;

import org.utgenome.util.TestHelper;
import org.utgenome.weaver.GenomeWeaver;
import org.xerial.util.FileUtil;
import org.xerial.util.log.Logger;

public class ConsensusBuilderTest
{
    private static Logger _logger = Logger.getLogger(ConsensusBuilderTest.class);

    @Test
    public void testname() throws Exception {
        File fasta = TestHelper.createTempFileFrom(ConsensusBuilderTest.class, "chr1-fragment.fa");
        File varDB = TestHelper.createTempFileFrom(ConsensusBuilderTest.class, "sample.freq");

        File newRef = FileUtil.createTempFile(new File("target"), "newref", ".fa");

        GenomeWeaver.execute(String.format("ConsensusBuilder -r %s -v %s -o %s", fasta, varDB, newRef));

        _logger.debug("done.");
    }
}
