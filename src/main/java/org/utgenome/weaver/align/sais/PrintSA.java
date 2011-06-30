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
// PrintSA.java
// Since: 2011/02/17
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.sais;

import java.io.StringWriter;

import org.utgenome.gwt.utgb.client.bio.IUPAC;
import org.utgenome.weaver.align.BWTransform;
import org.utgenome.weaver.align.GenomeWeaverCommand;
import org.utgenome.weaver.align.IUPACSequence;
import org.utgenome.weaver.align.LSeq;
import org.xerial.util.log.Logger;
import org.xerial.util.opt.Argument;
import org.xerial.util.opt.Option;

public class PrintSA extends GenomeWeaverCommand
{
    private static Logger _logger = Logger.getLogger(PrintSA.class);

    @Override
    public String name() {
        return "sa-debug";
    }

    @Override
    public String getOneLineDescription() {
        return "print SA and BWT for debugging";
    }

    @Override
    public Object getOptionHolder() {
        return this;
    }

    @Argument(index = 0, required = true)
    private String  seq;

    @Option(symbol = "r", description = "create SA for reverse string")
    private boolean isReverse     = false;

    @Option(symbol = "u", description = "Use Uint32SAIS")
    private boolean useUint32SAIS = false;

    @Override
    public void execute(String[] args) throws Exception {

        IUPACSequence s = new IUPACSequence(seq);
        if (isReverse)
            s = s.reverse();

        _logger.info("sequence: " + s);

        final int K = IUPAC.values().length;

        LSeq SA = new LSAIS.IntArray(new int[(int) s.textSize()], 0);
        LSAIS.suffixsort(s, SA, IUPAC.values().length);
        //LSeq SA = UInt32SAIS.SAIS(s, K);

        IUPACSequence bwt = new IUPACSequence(s.textSize());
        BWTransform.bwt(s, SA, bwt);

        //        WaveletArray wv = new WaveletArray(bwt, K);
        //        SparseSuffixArray ssa = SparseSuffixArray.buildFromSuffixArray(SA, 32);
        //        FMIndex fmIndex = new FMIndexOnWaveletArray(wv);

        for (int i = 0; i < SA.textSize(); ++i) {
            int sa = (int) SA.lookup(i);
            //            int sa_wv = (int) ssa.get(i, fmIndex);
            IUPAC c = bwt.getIUPAC(i);
            System.out.println(String.format("%3d %3d %s %s", i, sa, c == IUPAC.None ? "$" : c, rotateLeft(s, sa)));
        }
    }

    public static String rotateLeft(IUPACSequence s, int shift) {
        StringWriter w = new StringWriter();

        long len = s.textSize();
        for (int i = 0; i < len; ++i) {
            int index = i + shift;
            IUPAC c = s.getIUPAC(index % len);
            w.append(c.name());
        }

        return w.toString();
    }

}
