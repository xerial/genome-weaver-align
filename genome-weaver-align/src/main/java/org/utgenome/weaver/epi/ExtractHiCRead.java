/*------------------------CTGCTGTACCCTACATCCGCCTTGGCCGTACAGCAG--------------------------------------------------
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
// ExtractHiCRead.java
// Since: 2011/12/20
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.epi;

import java.io.BufferedInputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;

import org.utgenome.format.fastq.FastqRead;
import org.utgenome.format.fastq.FastqReader;
import org.utgenome.weaver.GenomeWeaverCommand;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.BitParallelSmithWaterman;
import org.utgenome.weaver.align.SmithWatermanAligner.Alignment;
import org.xerial.util.log.Logger;
import org.xerial.util.opt.Argument;
import org.xerial.util.opt.Option;

public class ExtractHiCRead extends GenomeWeaverCommand
{
    private static Logger _logger = Logger.getLogger(ExtractHiCRead.class);

    @Override
    public String getOneLineDescription() {
        return "Extract DNA sequence from HiC fastq data";
    }

    @Argument
    private List<File>   fastq                = new ArrayList<File>();

    @Option(symbol = "k", description = "num allowed mismatches")
    private int          numMismatchesAllowed = 3;

    @Option(symbol = "o")
    private String       outPrefix            = "fragment";

    private ACGTSequence adapter              = new ACGTSequence("CTGCTGTACCCTACATCCGCCTTGGCCGTACAGCAG");
    private ACGTSequence adapter_rc           = adapter.reverseComplement();

    @Override
    public void execute(String[] args) throws Exception {

        if (fastq.isEmpty()) {
            _logger.error("Two fastq files must be passed as the arguments");
            return;
        }

        InputStream in1 = new FileInputStream(fastq.get(0));
        if (fastq.get(0).getName().endsWith(".gz"))
            in1 = new GZIPInputStream(in1);

        InputStream in2 = new FileInputStream(fastq.get(1));
        if (fastq.get(1).getName().endsWith(".gz"))
            in2 = new GZIPInputStream(in2);

        FastqReader fin1 = new FastqReader(new InputStreamReader(new BufferedInputStream(in1)));
        FastqReader fin2 = new FastqReader(new InputStreamReader(new BufferedInputStream(in2)));

        Writer fout1 = new BufferedWriter(new FileWriter(String.format("%s_1.fastq", outPrefix)));
        Writer fout2 = new BufferedWriter(new FileWriter(String.format("%s_2.fastq", outPrefix)));

        try {

            while (true) {
                FastqRead f1 = fin1.next();
                FastqRead f2 = fin2.next();

                if (f1 == null || f2 == null) {
                    _logger.info("done.");
                    return;
                }

                ACGTSequence r1 = new ACGTSequence(f1.seq);
                ACGTSequence r2 = new ACGTSequence(f2.seq);

                Alignment adapterInF = findAdapter(r1);
                Alignment adapterInR = findAdapter(r2);

                if (adapterInF == null && adapterInR == null) {
                    _logger.warn("No common adapter found in %s and %s", f1.seqname, f2.seqname);
                    continue;
                }

                ACGTSequence dna1 = r1.subSequence(0, adapterInF.pos);
                ACGTSequence dna2 = r2.subSequence(0, adapterInR.pos);

                if (dna1.length() < 20 || dna2.length() < 20) {
                    _logger.warn("insufficient DNA length: %s, F:%d, R:%d : %s", f1.seqname, dna1.length(),
                            dna2.length(), f1.seq);
                    continue;
                }

                fout1.write(new FastqRead(f1.seqname, dna1.toString(), f1.qual.substring(0,
                        Math.min(adapterInF.pos, f1.qual.length()))).toFASTQString());
                fout2.write(new FastqRead(f2.seqname, dna2.toString(), f2.qual.substring(0,
                        Math.min(adapterInR.pos, f2.qual.length()))).toFASTQString());
            }
        }
        finally {
            fin1.close();
            fin2.close();

            fout1.close();
            fout2.close();
        }

    }

    private Alignment findAdapter(ACGTSequence ref) {
        Alignment aln = BitParallelSmithWaterman.alignBlockDetailed(ref, adapter, numMismatchesAllowed);
        if (aln == null) {
            aln = BitParallelSmithWaterman.alignBlockDetailed(ref, adapter_rc, numMismatchesAllowed);
        }
        return aln;
    }

}
