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
import java.io.FileOutputStream;
import java.io.IOException;

import org.utgenome.UTGBErrorCode;
import org.utgenome.UTGBException;
import org.utgenome.gwt.utgb.client.bio.IUPAC;
import org.utgenome.weaver.align.sais.LSAIS;
import org.utgenome.weaver.align.sais.UInt32Array;
import org.utgenome.weaver.align.sais.UInt32SAIS;
import org.xerial.util.StopWatch;
import org.xerial.util.log.Logger;
import org.xerial.util.opt.Argument;
import org.xerial.util.opt.Command;

/**
 * Performs burrows-wheeler transform
 * 
 * @author leo
 * 
 */
public class BWTransform implements Command
{
    private static Logger _logger = Logger.getLogger(BWTransform.class);

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

        StopWatch timer = new StopWatch();

        // Create IUPAC sequences (forward/reverse) from the given FASTA file
        Fasta2IUPAC.encode(fastaFile);

        BWTFiles forwardDB = new BWTFiles(fastaFile, Strand.FORWARD);
        BWTFiles reverseDB = new BWTFiles(fastaFile, Strand.REVERSE);

        // Create a suffix array and BWT string of the forward/reverse IUPAC sequence
        buildBWT(forwardDB);
        buildBWT(reverseDB);

        _logger.info(String.format("finised %.2f sec.", timer.getElapsedTime()));
    }

    public static void buildBWT(BWTFiles db) throws IOException, UTGBException {

        // Create BWT string
        pac2bwt(db);

        // Create a Wavelet array 
        {
            StopWatch timer = new StopWatch();
            _logger.info("Creating wavelet array " + db.bwtWavelet());
            IUPACSequence bwt = IUPACSequence.loadFrom(db.bwt());
            WaveletArray wv = new WaveletArray(bwt, 16);
            wv.saveTo(db.bwtWavelet());
            _logger.info(String.format("done. %.2f sec.", timer.getElapsedTime()));
        }

        //BWT2SparseSA.bwt2sparseSA(db);

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
                SA = new UInt32Array(seq.textSize());

            UInt32SAIS.SAIS(seq, SA, 16);
            //LSAIS.suffixsort(seq, SA, 16);
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

    public static void bwt(IUPACSequence seq, LSeq SA, IUPACSequence out) throws IOException {
        for (long i = 0; i < SA.textSize(); ++i) {
            if (SA.lookup(i) == 0) {
                out.setIUPAC(i, IUPAC.None);
            }
            else {
                out.setIUPAC(i, seq.getIUPAC(SA.lookup(i) - 1));
            }
        }
    }

}
