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

import org.utgenome.UTGBErrorCode;
import org.utgenome.UTGBException;
import org.xerial.util.log.Logger;
import org.xerial.util.opt.Argument;

public class Pac2BWT extends GenomeWeaverCommand
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

        BWTransform.pac2bwt(forwardDB);
        BWTransform.pac2bwt(reverseDB);
    }

}
