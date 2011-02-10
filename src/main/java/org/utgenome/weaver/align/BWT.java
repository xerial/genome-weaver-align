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
import java.io.File;
import java.io.FileOutputStream;

import org.utgenome.UTGBErrorCode;
import org.utgenome.UTGBException;
import org.utgenome.format.fasta.CompactFASTA;
import org.utgenome.format.fasta.FASTAPullParser;
import org.utgenome.gwt.utgb.client.bio.IUPAC;
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
        String outputFileName = FileType.removeFileExt(fastaFile) + ".iupac";
        _logger.info("input FASTA file: " + fastaFile);
        _logger.info("IUPAC file: " + outputFileName);

        BufferedOutputStream iupacFile = new BufferedOutputStream(new FileOutputStream(outputFileName));
        // Read the input FASTA file         
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
            _logger.info("sequence size: " + sequenceSize);
            offset = encoder.size();
        }
        encoder.close();
        _logger.info("total size: " + encoder.size());

    }

}
