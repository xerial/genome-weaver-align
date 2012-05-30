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
// AlignmentRecordTest.java
// Since: 2011/04/25
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.record;

import java.io.StringReader;

import org.junit.Test;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.record.AlignmentRecord;
import org.xerial.lens.SilkLens;
import org.xerial.util.log.Logger;

public class AlignmentRecordTest
{
    private static Logger _logger = Logger.getLogger(AlignmentRecordTest.class);

    @Test
    public void serialize() throws Exception {

        AlignmentRecord r = new AlignmentRecord();
        r.chr = "chr1";
        r.start = 43;
        r.strand = Strand.FORWARD;
        r.numMismatches = 0;
        r.setCIGAR("10S90M");

        String silk = SilkLens.toSilk(r);
        _logger.info(silk);
        AlignmentRecord r2 = SilkLens.loadSilk(AlignmentRecord.class, new StringReader(silk));

        _logger.info(SilkLens.toSilk(r2));

    }
}
