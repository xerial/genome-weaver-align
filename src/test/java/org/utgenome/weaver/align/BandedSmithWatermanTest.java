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

import org.junit.Test;
import org.utgenome.weaver.align.SmithWatermanAligner.Alignment;
import org.xerial.util.log.Logger;

public class BandedSmithWatermanTest
{
    private static Logger _logger = Logger.getLogger(BandedSmithWatermanTest.class);

    public void align(String ref, String query) throws Exception {
        ACGTSequence r = new ACGTSequence(ref);
        ACGTSequence q = new ACGTSequence(query);
        Alignment alignment = BandedSmithWaterman.align(r, q);
        _logger.debug(alignment);
    }

    public void align(String ref, String query, AlignmentScoreConfig config) throws Exception {
        ACGTSequence r = new ACGTSequence(ref);
        ACGTSequence q = new ACGTSequence(query);
        Alignment alignment = BandedSmithWaterman.align(r, q, config);
        _logger.debug(alignment);
    }

    @Test
    public void align() throws Exception {
        align("GATATAGAGATCTGGCCTAG", "TAGAGAT");
    }

    @Test
    public void alignGap() throws Exception {
        AlignmentScoreConfig config = new AlignmentScoreConfig();
        config.gapOpenPenalty = 3;
        config.gapExtensionPenalty = 5;
        align("TATACCAAGATATAGAGATCTGGCAAGTGTGTTAT", "TATACCAAGATATAGATCTGGCAAGTGTGTTAT", config);
    }
}
