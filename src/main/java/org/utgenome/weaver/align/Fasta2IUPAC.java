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
// EncodeFasta.java
// Since: 2011/02/17
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;

import org.utgenome.format.fasta.CompactFASTA;
import org.utgenome.format.fasta.FASTAPullParser;
import org.utgenome.gwt.utgb.client.bio.IUPAC;
import org.utgenome.weaver.GenomeWeaverCommand;
import org.utgenome.weaver.align.SequenceBoundary.SequenceIndex;
import org.xerial.lens.SilkLens;
import org.xerial.silk.SilkWriter;
import org.xerial.util.log.Logger;
import org.xerial.util.opt.Argument;

public class Fasta2IUPAC extends GenomeWeaverCommand
{
    private static Logger _logger = Logger.getLogger(Fasta2IUPAC.class);

    @Override
    public String name() {
        return "fasta2iupac";
    }

    @Override
    public String getOneLineDescription() {
        return "Fasta -> iupac binary";
    }

    @Override
    public Object getOptionHolder() {
        return this;
    }

    @Argument(index = 0, required = true)
    private String fastaFile;

    @Override
    public void execute(String[] args) throws Exception {

        encode(fastaFile);
    }

    public static void encode(String fastaFile) throws IOException {
        encodeFASTA(fastaFile);

        BWTFiles forwardDB = new BWTFiles(fastaFile, Strand.FORWARD);
        BWTFiles reverseDB = new BWTFiles(fastaFile, Strand.REVERSE);
        reverse(forwardDB, reverseDB);
    }

    protected static void encodeFASTA(String fastaFile) throws IOException {
        BWTFiles forwardDB = new BWTFiles(fastaFile, Strand.FORWARD);
        _logger.info("input FASTA file: " + fastaFile);

        long totalSize = -1;
        {
            // Read the input FASTA file, then encode the sequences using the IUPAC code
            BufferedOutputStream iupacOut = new BufferedOutputStream(new FileOutputStream(forwardDB.iupac()));
            SilkWriter indexOut = new SilkWriter(new BufferedWriter(new FileWriter(forwardDB.pacIndex())));
            IUPACSequenceWriter encoder = new IUPACSequenceWriter(iupacOut);
            FASTAPullParser fasta = new FASTAPullParser(new File(fastaFile));
            long lineCount = 1;
            long offset = 0;
            for (String desc; (desc = fasta.nextDescriptionLine()) != null; lineCount++) {
                String seqName = CompactFASTA.pickSequenceName(desc);
                _logger.info(String.format("reading %s", seqName));
                for (String seq; (seq = fasta.nextSequenceLine()) != null; lineCount++) {
                    seq = seq.trim();
                    for (int i = 0; i < seq.length(); ++i) {
                        // 'A' .. 'Z'
                        char base = Character.toUpperCase(seq.charAt(i));
                        IUPAC iupac = IUPAC.encode(base);
                        if (iupac == IUPAC.None) {
                            // illegal character
                            _logger.warn(String.format("illegal character '%s' at line:%,d, pos:%d: Use N instead",
                                    base, lineCount, i + 1));
                            iupac = IUPAC.N;
                        }
                        encoder.append(iupac);
                    }
                }
                long pos = encoder.size();
                long sequenceSize = pos - offset;
                SequenceIndex index = new SequenceIndex(seqName, desc, sequenceSize, offset);
                indexOut.leafObject("index", index);
                _logger.debug("\n" + SilkLens.toSilk("index", index));
                offset = encoder.size();
            }
            // This part ensure the sequence size (including sentinel $) becomes 2n, which fits in a byte array 
            if (encoder.size() % 2 == 0)
                encoder.append(IUPAC.N);
            // append a sentinel
            encoder.append(IUPAC.None);

            encoder.close();
            totalSize = encoder.size();
            _logger.info(String.format("total num bases: %,d", totalSize));
            indexOut.leaf("total size", totalSize);
            indexOut.close();
        }

    }

    protected static void reverse(BWTFiles forwardDB, BWTFiles reverseDB) throws IOException {
        // Reverse the IUPAC sequence
        IUPACSequence forwardSeq = IUPACSequence.loadFrom(forwardDB.iupac());
        _logger.info("Reverse the sequence");
        _logger.info("Reverse IUPAC file: " + reverseDB.iupac());
        IUPACSequenceWriter encoder = new IUPACSequenceWriter(new BufferedOutputStream(new FileOutputStream(
                reverseDB.iupac())));
        // reverse IN[0..n-2] (excludes the sentinel)
        long i = forwardSeq.textSize() - 2;
        //        if (forwardSeq.getIUPAC(i) == IUPAC.N) {
        //            i--;
        //        }
        for (; i >= 0; --i) {
            encoder.append(forwardSeq.getIUPAC(i));
        }
        // append a sentinel.
        encoder.append(IUPAC.None);
        encoder.close();
    }

}
