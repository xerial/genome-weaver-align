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

import java.io.IOException;

import org.utgenome.UTGBErrorCode;
import org.utgenome.UTGBException;
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
        EncodeFasta.encode(fastaFile);

        BWTFiles forwardDB = new BWTFiles(fastaFile, Strand.FORWARD);
        BWTFiles reverseDB = new BWTFiles(fastaFile, Strand.REVERSE);

        // Create a suffix array and BWT string of the forward/reverse IUPAC sequence
        buildBWT(forwardDB);
        buildBWT(reverseDB);

        _logger.info(String.format("finised %.2f sec.", timer.getElapsedTime()));
    }

    public static void buildBWT(BWTFiles db) throws IOException, UTGBException {

        // Create BWT string
        Pac2BWT.pac2bwt(db);

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

}
