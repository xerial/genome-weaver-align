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
// BWT.java
// Since: 2011/02/10
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
import java.util.Arrays;

import org.utgenome.UTGBErrorCode;
import org.utgenome.UTGBException;
import org.utgenome.format.fasta.CompactFASTA;
import org.utgenome.format.fasta.FASTAPullParser;
import org.utgenome.gwt.utgb.client.bio.IUPAC;
import org.utgenome.weaver.align.IUPACSequence.SequenceIndex;
import org.xerial.lens.SilkLens;
import org.xerial.silk.SilkWriter;
import org.xerial.util.FileType;
import org.xerial.util.log.Logger;
import org.xerial.util.opt.Argument;
import org.xerial.util.opt.Command;

/**
 * Performs burrows-wheeler transform
 * 
 * @author leo
 * 
 */
public class BWT implements Command
{
    private static Logger _logger = Logger.getLogger(BWT.class);

    @Override
    public String name() {
        return "bwt";
    }

    @Override
    public String getOneLineDescription() {
        return "Burrows-Wheeler Transformation (BWT)";
    }

    @Override
    public Object getOptionHolder() {
        return this;
    }

    /**
     * input FASTA file (.fa, .tar.gz, .fa.gz, types are allowed)
     */
    @Argument(index = 0)
    private String fastaFile;

    @Override
    public void execute(String[] args) throws Exception {

        if (fastaFile == null)
            throw new UTGBException(UTGBErrorCode.MISSING_FILES, "no input FASTA file is given");

        // Output IUPAC sequence to a file
        final String fastaPrefix = FileType.removeFileExt(fastaFile);
        String iupacFileName = fastaPrefix + ".iupac";
        String indexFileName = fastaPrefix + ".i.silk";
        _logger.info("input FASTA file: " + fastaFile);
        _logger.info("IUPAC file: " + iupacFileName);
        _logger.info("index file: " + indexFileName);

        {
            BufferedOutputStream iupacFile = new BufferedOutputStream(new FileOutputStream(iupacFileName));
            SilkWriter indexOut = new SilkWriter(new BufferedWriter(new FileWriter(indexFileName)));
            // Read the input FASTA file, then encode the sequences using the IUPAC code         
            IUPACSequenceWriter encoder = new IUPACSequenceWriter(iupacFile);
            FASTAPullParser fasta = new FASTAPullParser(new File(fastaFile));
            int lineCount = 1;
            int offset = 0;
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
                            _logger.warn(String.format("illegal character '%s' at line:%,d, char:%d, char:%s", base,
                                    lineCount, i + 1));
                            continue;
                        }

                        encoder.append(iupac);
                    }
                }

                int pos = encoder.size();
                int sequenceSize = pos - offset;
                SequenceIndex index = new SequenceIndex(seqName, desc, sequenceSize, offset);
                indexOut.leafObject("index", index);
                _logger.info("\n" + SilkLens.toSilk("index", index));
                offset = encoder.size();
            }
            encoder.close();
            _logger.info("total size: " + encoder.size());
            indexOut.leaf("total size", encoder.size());

            indexOut.close();
        }

        // Reverse the IUPAC sequence
        IUPACSequence forwardSeq = new IUPACSequence(new File(iupacFileName));
        {
            String reverseIupacFileName = fastaPrefix + ".r.iupac";
            _logger.info("Reverse the sequence");
            _logger.info("Reverse IUPAC file: " + reverseIupacFileName);

            BufferedOutputStream revOut = new BufferedOutputStream(new FileOutputStream(reverseIupacFileName));
            forwardSeq.reverse(revOut);
            revOut.close();
        }

        {
            // Create a suffix array of the forward IUPAC sequence
            int[] SA = new int[forwardSeq.size()];
            SAIS.suffixsort(forwardSeq, SA, 16);
            File suffixArrayFile = new File(fastaPrefix + ".sa");
            _logger.info("SA file: " + suffixArrayFile);
            BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(suffixArrayFile));
            final long byteSize = SA.length * 4;
            for (int i = 0; i < SA.length; ++i) {
                out.write((SA[i] >>> 24) & 0xFF);
                out.write((SA[i] >>> 16) & 0xFF);
                out.write((SA[i] >>> 8) & 0xFF);
                out.write((SA[i]) & 0xFF);
            }
            out.close();

            // Create a BWT string of the forward IUPAC sequence from the generated suffix array
            IUPAC[] bwt = bwt(forwardSeq, SA);
            _logger.info("SA : " + Arrays.toString(SA));
            _logger.info("IN : " + Arrays.toString(forwardSeq.toArray()));
            _logger.info("BWT: " + Arrays.toString(bwt));
        }

    }

    public static IUPAC[] bwt(IUPACSequence seq, int[] SA) {
        IUPAC[] bwt = new IUPAC[SA.length];

        for (int i = 0; i < SA.length; ++i) {
            if (SA[i] == 0) {
                bwt[i] = IUPAC.None;
            }
            else {
                bwt[i] = seq.getIUPAC(SA[i] - 1);
            }
        }

        return bwt;
    }

}
