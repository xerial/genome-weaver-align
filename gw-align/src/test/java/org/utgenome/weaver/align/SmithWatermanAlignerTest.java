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
// BandedSmithWatermanTest.java
// Since: 2011/08/29
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import static org.junit.Assert.*;

import org.junit.Test;
import org.utgenome.weaver.align.SmithWatermanAligner.Alignment;
import org.xerial.util.log.Logger;

public class SmithWatermanAlignerTest
{
    private static Logger _logger = Logger.getLogger(SmithWatermanAlignerTest.class);

    public Alignment bandedAlign(String ref, String query) throws Exception {
        ACGTSequence r = new ACGTSequence(ref);
        ACGTSequence q = new ACGTSequence(query);
        Alignment alignment = SmithWatermanAligner.bandedAlign(r, q);
        _logger.debug(alignment);
        return alignment;
    }

    public Alignment bandedAlign(String ref, String query, AlignmentScoreConfig config) throws Exception {
        ACGTSequence r = new ACGTSequence(ref);
        ACGTSequence q = new ACGTSequence(query);
        Alignment alignment = SmithWatermanAligner.bandedAlign(r, q, config);
        _logger.debug(alignment);
        return alignment;
    }

    @Test
    public void align() throws Exception {
        Alignment alignment = bandedAlign("GATATAGAGATCTGGCCTAG", "TAGAGAT");
        assertEquals("7M", alignment.cigar.toString());
        assertEquals(4, alignment.pos);
        assertEquals(0, alignment.numMismatches);
    }

    @Test
    public void alignMismatch() throws Exception {
        Alignment alignment = bandedAlign("GATATAGAGATCTGGCCTAG", "TATACAGATCTGGCCTAG");
        assertEquals("4M1X13M", alignment.cigar.toString());
        assertEquals(2, alignment.pos);
        assertEquals(1, alignment.numMismatches);
    }

    @Test
    public void alignInsertion() throws Exception {
        AlignmentScoreConfig config = new AlignmentScoreConfig();
        config.gapOpenPenalty = 5;
        config.gapExtensionPenalty = 2;
        Alignment alignment = bandedAlign("TATACCAAGATATAGATCTGGCAAGTGTGTTAT", "TACCAAGATATAGAGATCTGGCAAGTGTGTTAT",
                config);
        assertEquals("11M2I20M", alignment.cigar.toString());
        assertEquals(2, alignment.pos);
        assertEquals(2, alignment.numMismatches);
    }

    @Test
    public void alignDeletion() throws Exception {
        AlignmentScoreConfig config = new AlignmentScoreConfig();
        config.gapOpenPenalty = 5;
        config.gapExtensionPenalty = 2;
        Alignment alignment = bandedAlign("TATACCAAGATATAGAGATCTGGCAAGTGTGTTAT", "TACCAAGATATAGATCTGGCAAGTGTGTTAT",
                config);
        assertEquals("11M2D20M", alignment.cigar.toString());
        assertEquals(2, alignment.pos);
        assertEquals(2, alignment.numMismatches);
    }

    @Test
    public void clippedAlignment() throws Exception {
        Alignment alignment = bandedAlign("ATTTGTATTATACCAAGAT", "ACCGGAAGGACCAAGAT");
        assertEquals("9S8M", alignment.cigar.toString());
        assertEquals(11, alignment.pos);
        assertEquals(0, alignment.numMismatches);
    }
}
