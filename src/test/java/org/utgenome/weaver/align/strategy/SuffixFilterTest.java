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
// SuffixFilterTest.java
// Since: 2011/07/27
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.strategy;

import org.junit.Test;
import org.utgenome.UTGBException;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.SequenceBoundary.PosOnGenome;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.strategy.SuffixFilter.Candidate;
import org.utgenome.weaver.parallel.Reporter;
import org.xerial.lens.SilkLens;
import org.xerial.util.log.Logger;

public class SuffixFilterTest
{
    private static Logger _logger = Logger.getLogger(SuffixFilterTest.class);

    @Test
    public void constructor() throws Exception {
        final FMIndexOnGenome fmIndex = FMIndexOnGenome.buildFromSequence("seq",
                "TAGCCTATAGAGCGAAAGGAGATATATAGCCCGAGTAT");

        final ACGTSequence q = new ACGTSequence("GCCTATAGAGCG");
        SuffixFilter f = new SuffixFilter(2, fmIndex, q, Strand.FORWARD);
        f.match(new Reporter() {
            @Override
            public void emit(Object result) throws UTGBException {
                Candidate input = (Candidate) result;
                _logger.debug(SilkLens.toSilk("match", input));
                PosOnGenome gc = fmIndex.toGenomeCoordinate(input.si.lowerBound, input.offset, Strand.FORWARD);
                if (gc != null)
                    _logger.debug(SilkLens.toSilk("loc", gc));
            }

        });
    }

    @Test
    public void align() throws Exception {
        final FMIndexOnGenome fmIndex = FMIndexOnGenome.buildFromSequence("seq", "ATATAGCCCGAGTAT");

        final ACGTSequence q = new ACGTSequence("TATAGCCC");
        SuffixFilter f = new SuffixFilter(2, fmIndex, q, Strand.FORWARD);
        f.match(new Reporter() {
            @Override
            public void emit(Object result) throws UTGBException {
                Candidate input = (Candidate) result;
                _logger.debug(SilkLens.toSilk("match", input));
                PosOnGenome gc = fmIndex.toGenomeCoordinate(input.si.lowerBound, input.offset, Strand.FORWARD);
                if (gc != null)
                    _logger.debug(SilkLens.toSilk("loc", gc));
            }

        });

    }

}
