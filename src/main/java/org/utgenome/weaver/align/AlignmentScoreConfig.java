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

/**
 * Alignment scoring configuration
 * 
 * @author leo
 * 
 */
public class AlignmentScoreConfig
{
    public final int numMismatchesAllowed = 2;

    public final int matchScore           = 1;
    public final int mismatchPenalty      = 3;
    public final int gapOpenPenalty       = 11;
    public final int gapExtentionPenalty  = 4;

    public final int numGapOpenAllowed    = 1;
    public final int numSplitAlowed       = 1;

}
