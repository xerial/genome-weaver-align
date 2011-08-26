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
// AlignmentConfig.java
// Since: 2011/08/16
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.util.List;

import org.xerial.util.opt.Argument;
import org.xerial.util.opt.Option;

/**
 * Alignment command configuration
 * 
 * @author leo
 * 
 */
public class AlignmentConfig extends AlignmentScoreConfig
{
    @Option(symbol = "r", description = "reference sequence")
    public String       refSeq;

    @Option(symbol = "q", description = "single query sequence")
    public String       query;

    @Argument(name = "read file")
    public List<String> readFiles;

    public static enum Strategy {
        SF("suffix filter"), BD("bi-directional search"), BWA("bwa");
        public final String description;

        private Strategy(String description) {
            this.description = description;
        }
    }

    @Option(symbol = "m", description = "alignment strategy. sf(suffix filter), bd(bidirectinal search), bwa (best-hit first)")
    public Strategy strategy = Strategy.SF;

    public static enum ReportType {
        BESTHIT, ALLHITS, TOPL
    }

    @Option(symbol = "R", description = "reporting method. besthit (default), allhits, topL (top-L hits)")
    public ReportType reportType;

    @Option(symbol = "L", description = "number of hits to report (default=5). Used only when -R topL option is set")
    public int        topL = 5;

}
