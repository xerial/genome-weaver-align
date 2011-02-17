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
// IUPAC2BWT.java
// Since: 2011/02/16
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;

import org.utgenome.UTGBErrorCode;
import org.utgenome.UTGBException;
import org.utgenome.gwt.utgb.client.bio.IUPAC;
import org.xerial.util.StopWatch;
import org.xerial.util.log.Logger;
import org.xerial.util.opt.Argument;
import org.xerial.util.opt.Command;

public class Pac2BWT implements Command
{
    private static Logger _logger = Logger.getLogger(Pac2BWT.class);

    @Override
    public String name() {
        return "pac2bwt";
    }

    @Override
    public String getOneLineDescription() {
        return "create BWT from an IUPAC file";
    }

    @Override
    public Object getOptionHolder() {
        return this;
    }

    @Argument(index = 0)
    private String fastaFile;

    @Override
    public void execute(String[] args) throws Exception {
        if (fastaFile == null)
            throw new UTGBException(UTGBErrorCode.MISSING_FILES, "no input fasta file");

        BWTFiles forwardDB = new BWTFiles(fastaFile, Strand.FORWARD);
        BWTFiles reverseDB = new BWTFiles(fastaFile, Strand.REVERSE);

        pac2bwt(forwardDB);
        pac2bwt(reverseDB);
    }

    public static void pac2bwt(BWTFiles db) throws UTGBException, IOException {

        StopWatch timer = new StopWatch();
        {
            IUPACSequence seq = IUPACSequence.loadFrom(db.iupac());
            timer.reset();
            _logger.info("Creating a suffix array of " + db.iupac());
            LSeq SA = null;
            if (seq.textSize() < Integer.MAX_VALUE) {
                SA = new LSAIS.IntArray(new int[(int) seq.textSize()], 0);
            }
            else
                SA = new LIntArray(seq.textSize());

            LSAIS.suffixsort(seq, SA, 16);
            _logger.info(String.format("%.2f sec.", timer.getElapsedTime()));

            _logger.info("Creating sparse suffix array " + db.sparseSuffixArray());
            SparseSuffixArray ssa = SparseSuffixArray.buildFromSuffixArray(SA, 32);
            ssa.saveTo(db.sparseSuffixArray());

            _logger.info("Creating a BWT string: " + db.bwt());
            timer.reset();
            IUPACSequenceWriter writer = new IUPACSequenceWriter(new BufferedOutputStream(
                    new FileOutputStream(db.bwt())));
            bwt(seq, SA, writer);
            writer.close();
            _logger.info(String.format("%.2f sec.", timer.getElapsedTime()));
        }

    }

    public static void bwt(IUPACSequence seq, LSeq SA, IUPACSequenceWriter out) throws IOException {
        for (long i = 0; i < SA.textSize(); ++i) {
            if (SA.lookup(i) == 0) {
                out.append(IUPAC.None);
            }
            else {
                out.append(seq.getIUPAC(SA.lookup(i) - 1));
            }
        }
    }

}
