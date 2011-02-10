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
// BWTTest.java
// Since: 2011/02/10
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.io.File;

import org.junit.Test;
import org.utgenome.util.TestHelper;
import org.utgenome.weaver.GenomeWeaver;

public class BWTTest
{
    @Test
    public void bwt() throws Exception {
        File fastaArchive = TestHelper.createTempFileFrom(BWTTest.class, "sample-archive.fa.tar.gz");
        GenomeWeaver.execute(String.format("bwt %s", fastaArchive));
    }
}
