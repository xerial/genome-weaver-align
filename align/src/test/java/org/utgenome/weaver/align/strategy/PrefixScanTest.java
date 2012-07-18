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
// PrefixScanTest.java
// Since: 2011/09/15
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.strategy;

import org.junit.Before;
import org.junit.Test;
import org.utgenome.format.fasta.FASTAPullParser;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.PackFasta;
import org.utgenome.weaver.align.PackFasta.PackedFASTA;
import org.utgenome.weaver.align.PackFastaTest;
import org.utgenome.weaver.align.Strand;
import org.xerial.util.FileResource;
import org.xerial.util.log.Logger;

public class PrefixScanTest
{
    private static Logger   _logger = Logger.getLogger(PrefixScanTest.class);

    private FMIndexOnGenome fmIndex;

    @Before
    public void setUp() throws Exception {
        PackedFASTA fasta = PackFasta.encode(new FASTAPullParser(FileResource.open(PackFastaTest.class, "sample2.fa")));
        fmIndex = FMIndexOnGenome.buildFromSequence("sample", fasta.sequence);
    }

    @Test
    public void scan() throws Exception {
        ACGTSequence q = new ACGTSequence(
                "TTTATACGTGTTTGAATAACTTGGCCAAATCGCCGAGAAGGAATAGAATACTGGACGACATTGTACATATTTTCCAAAAAATCAGAAAGTAGATGACGGG");
        PrefixScan scan = PrefixScan.scanRead(fmIndex, q, Strand.FORWARD, new StaircaseFilter(q.length(), 7));
        _logger.debug(scan);
    }
}
