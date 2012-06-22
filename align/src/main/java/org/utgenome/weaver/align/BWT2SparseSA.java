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
// BWT2SparseSA.java
// Since: 2011/02/17
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.io.IOException;

import org.utgenome.UTGBErrorCode;
import org.utgenome.UTGBException;
import org.utgenome.weaver.GenomeWeaverCommand;
import org.xerial.util.StopWatch;
import org.xerial.util.log.Logger;
import org.xerial.util.opt.Argument;

public class BWT2SparseSA extends GenomeWeaverCommand
{
    private static Logger _logger = Logger.getLogger(BWT2SparseSA.class);

    @Override
    public String name() {
        return "bwt2sa";
    }

    @Override
    public String getOneLineDescription() {
        return "Create a sparse suffix array from BWT string";
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
            throw new UTGBException(UTGBErrorCode.MISSING_OPTION, "no FASTA file is given");

        BWTFiles forwardDB = new BWTFiles(fastaFile, Strand.FORWARD);
        BWTFiles reverseDB = new BWTFiles(fastaFile, Strand.REVERSE);
        bwt2sparseSA(forwardDB);
        bwt2sparseSA(reverseDB);
    }

    public static void bwt2sparseSA(BWTFiles db) throws IOException {
        StopWatch timer = new StopWatch();
        {
            _logger.info("Loading bwt string: " + db.bwt());
            ACGTSequence bwt = ACGTSequence.loadFrom(db.bwt());
            _logger.info("Creating sparse suffix array " + db.sparseSuffixArray());
            FMIndex fm = new FMIndexOnOccTable(bwt);
            _logger.info("Done.");
            SparseSuffixArray ssa = SparseSuffixArray.createFromBWT(fm, 32);
            ssa.saveTo(db.sparseSuffixArray());
        }
        _logger.info(String.format("%.2f sec.", timer.getElapsedTime()));
    }

}
