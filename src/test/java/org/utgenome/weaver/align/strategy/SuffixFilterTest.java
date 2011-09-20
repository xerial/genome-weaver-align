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
// SuffixFilterTest.java
// Since: 2011/09/15
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.strategy;

import java.util.List;

import org.junit.BeforeClass;
import org.junit.Test;
import org.utgenome.format.fasta.FASTAPullParser;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.AlignmentConfig;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.PackFasta;
import org.utgenome.weaver.align.PackFasta.PackedFASTA;
import org.utgenome.weaver.align.PackFastaTest;
import org.utgenome.weaver.align.record.AlignmentRecord;
import org.xerial.util.FileResource;

public class SuffixFilterTest
{
    private static ACGTSequence    ref;
    private static FMIndexOnGenome fmIndex;
    private static AlignmentConfig config = new AlignmentConfig();

    @BeforeClass
    public static void setUp() throws Exception {
        PackedFASTA fasta = PackFasta.encode(new FASTAPullParser(FileResource.open(PackFastaTest.class, "sample3.fa")));
        ref = fasta.sequence;
        fmIndex = FMIndexOnGenome.buildFromSequence("sample", fasta.sequence);
    }

    public static AlignmentRecord align(String query) throws Exception {
        return align(new ACGTSequence(query));
    }

    public static AlignmentRecord align(ACGTSequence q) throws Exception {
        SuffixFilter f = new SuffixFilter(fmIndex, ref, config);
        List<AlignmentRecord> result = f.align(q);
        if (result.size() == 0)
            return null;
        else
            return result.get(0);
    }

    @Test
    public void exactMatch() throws Exception {
        align("CACTTTAGTATAATTGTTTTTAGTTTTTGGCAAAACTATTGTCTAAACAG");

    }

    @Test
    public void insertion() throws Exception {
        align("CACTTTAGTATAATTGTTTTTAGCCTTTTTGGCAAAACTATTGTCTAAACAG");

    }

}
