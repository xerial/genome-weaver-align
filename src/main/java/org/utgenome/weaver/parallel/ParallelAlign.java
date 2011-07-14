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
// ParallelAlign.java
// Since: 2011/05/25
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.parallel;

import org.utgenome.weaver.align.AlignmentScoreConfig;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.record.AlignmentRecord;
import org.utgenome.weaver.align.record.ReadSequence;
import org.utgenome.weaver.align.record.ReadSequenceReader;
import org.utgenome.weaver.align.record.ReadSequenceReaderFactory;

public class ParallelAlign
{
    private FMIndexOnGenome      fmIndex;
    private AlignmentScoreConfig config;

    public ParallelAlign(FMIndexOnGenome fmIndex, AlignmentScoreConfig config) {
        this.fmIndex = fmIndex;
        this.config = config;
    }

    public void map(Provider<ReadSequence> readSet, Reporter reporter) {
        ReadSequence seq = readSet.get();
        // do some alignment
        AlignmentRecord rec = new AlignmentRecord();
        reporter.emit(rec);
    }

    public static void execute(String fastaFilePrefix, String fastqFile) throws Exception {

        FMIndexOnGenome fmIndex = new FMIndexOnGenome(fastaFilePrefix);
        AlignmentScoreConfig config = new AlignmentScoreConfig();

        ParallelAlign cmd = new ParallelAlign(fmIndex, config);
        ReadSequenceReader fastq = ReadSequenceReaderFactory.createFASTQReader(fastqFile);

    }
}
