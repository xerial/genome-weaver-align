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
// AlignmentScoreConfig.java
// Since: 2011/04/28
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import org.xerial.util.opt.Option;

/**
 * Alignment scoring configuration
 * 
 * @author leo
 * 
 */
public class AlignmentScoreConfig
{
    @Option(symbol = "k", description = "maximum edit distances")
    public int maximumEditDistances   = 2;
    @Option(symbol = "g", description = "# of gap open allowed. default=1")
    public int numGapOpenAllowed      = 1;
    @Option(symbol = "e", description = "# of gap extension allowed. default=4")
    public int numGapExtensionAllowed = 4;
    @Option(symbol = "s", description = "# of read split allowed. default=1")
    public int numSplitAlowed         = 1;

    @Option(symbol = "M", description = "match score. default=1")
    public int matchScore             = 1;
    @Option(symbol = "N", description = "mismatch penalty. default=3")
    public int mismatchPenalty        = 3;
    @Option(symbol = "G", description = "gap open penalty. default=11")
    public int gapOpenPenalty         = 11;
    @Option(symbol = "E", description = "gap extension penalty. default=4")
    public int gapExtensionPenalty    = 4;
    @Option(symbol = "S", description = "split open penalty. default=5")
    public int splitOpenPenalty       = 5;

    @Option(symbol = "P", description = "skip indels within P bases from read head and tail")
    public int indelEndSkip           = 5;

}
