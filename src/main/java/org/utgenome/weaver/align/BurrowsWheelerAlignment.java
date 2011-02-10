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
// BurrowsWheelerAlignment.java
// Since: 2011/02/10
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.io.File;

import org.utgenome.gwt.utgb.client.bio.IUPAC;
import org.utgenome.weaver.align.IUPACSequence.IUPACBinaryInfo;
import org.xerial.util.opt.Argument;
import org.xerial.util.opt.Command;
import org.xerial.util.opt.Option;

/**
 * Burrows-Wheeler aligner for IUPAC sequences
 * 
 * @author leo
 * 
 */
public class BurrowsWheelerAlignment implements Command
{

    @Override
    public String name() {
        return "align";
    }

    @Override
    public String getOneLineDescription() {
        return "performs alignment";
    }

    @Override
    public Object getOptionHolder() {
        return this;
    }

    @Argument(index = 0)
    private String fastaFilePrefix;

    @Option(symbol = "L", description = "Save 1/L for occurrence count table (default = 256)")
    private int    L = 256;

    @Override
    public void execute(String[] args) throws Exception {

        IUPACBinaryInfo index = IUPACBinaryInfo.loadSilk(IUPACBinaryInfo.getFileName(fastaFilePrefix));
        final int N = index.totalSize;
        final int K = IUPAC.values().length;

        // loading BWT sequences
        File bwtForwardFile = new File(fastaFilePrefix + ".bwt");
        File bwtReverseFile = new File(fastaFilePrefix + ".rbwt");
        IUPACSequence bwtF = new IUPACSequence(bwtForwardFile, N);
        IUPACSequence bwtR = new IUPACSequence(bwtReverseFile, N);

        // Compute the occurrence tables
        OccurrenceCountTable occF = new OccurrenceCountTable(bwtF, L);
        OccurrenceCountTable occR = new OccurrenceCountTable(bwtR, L);

        // Compute the character start index
        CharacterCount C = new CharacterCount(bwtF);

    }
}
