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
// SmithWatermanAlign.java
// Since: 2011/08/30
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.utgenome.UTGBException;
import org.utgenome.gwt.utgb.client.bio.ChrLoc;
import org.utgenome.shell.UTGBShellException;
import org.utgenome.weaver.GenomeWeaverCommand;
import org.utgenome.weaver.align.SmithWatermanAligner.Alignment;
import org.xerial.util.log.Logger;
import org.xerial.util.opt.Argument;
import org.xerial.util.opt.Option;

public class SmithWatermanAlign extends GenomeWeaverCommand
{
    private static Logger _logger = Logger.getLogger(SmithWatermanAlign.class);

    @Override
    public String name() {
        return "sw-align";
    }

    @Override
    public String getOneLineDescription() {
        return "Simple Smith Waterman alignment";
    }

    @Override
    public Object getOptionHolder() {
        return config;
    }

    public static class Config extends AlignmentScoreConfig
    {

        @Option(symbol = "r", description = "reference name")
        public String reference;

        @Argument(index = 0, name = "query range (e.g., chr1:100-200)")
        public String range;

        @Option(symbol = "q", description = "query sequence")
        public String query;

    }

    public static class RegionQueryExpr
    {
        private static Pattern p = Pattern.compile("([^:]+)(:([0-9]+)(-([0-9]+))?)?");

        public static ChrLoc parse(String expr) throws UTGBShellException {
            Matcher m = p.matcher(expr.replaceAll(",", "")); // remove comma
            if (!m.matches())
                throw new UTGBShellException("invalid query format:" + expr);
            String chr = m.group(1);
            String sStart = m.group(3);
            String sEnd = m.group(5);

            int start = 0;
            if (sStart != null)
                start = Integer.parseInt(sStart);
            int end = Integer.MAX_VALUE;
            if (sEnd != null)
                end = Integer.parseInt(sEnd);
            else
                end = start;

            return new ChrLoc(chr, start, end);
        }
    }

    private Config config = new Config();

    @Override
    public void execute(String[] args) throws Exception {

        if (config.reference == null)
            throw new UTGBException("No reference is given. Use -r option");

        if (config.query == null)
            throw new UTGBException("No query is given. Use -q option");

        if (config.range == null)
            throw new UTGBException("No range is given. Specify in chr:start-end format. (e.g. chr1:1000-2000)");

        ChrLoc loc = RegionQueryExpr.parse(config.range);

        _logger.info("Loading reference: %s", config.reference);
        BWTFiles bwtFiles = new BWTFiles(config.reference, Strand.FORWARD);
        ACGTSequence ref = ACGTSequence.loadFrom(bwtFiles.pac());
        SequenceBoundary boundary = SequenceBoundary.loadSilk(bwtFiles.pacIndex());

        ACGTSequence query = new ACGTSequence(config.query);
        int range = Math.max(loc.length(), query.length());
        long offset = boundary.toIndex(loc.chr, loc.start);
        ACGTSequence target = ref.subSequence(offset, offset + range);
        Alignment alignF = SmithWatermanAligner.align(target, query, config, false);
        Alignment alignR = SmithWatermanAligner.align(target, query.reverseComplement(), config, false);
        System.out.println(String.format("query: %s:%,d-%,d", loc.chr, loc.start, loc.start + range));

        if (alignF == null && alignR == null)
            System.out.println("No match");
        else {
            if (alignF.numMismatches <= alignR.numMismatches) {
                reportAlignment(loc, alignF, Strand.FORWARD);
            }
            else {
                reportAlignment(loc, alignR, Strand.REVERSE);
            }
        }
    }

    private void reportAlignment(ChrLoc loc, Alignment aln, Strand strand) {
        System.out.println(String.format("Match at %s:%d %s", loc.chr, loc.start + aln.pos, strand));
        System.out.println(aln.toString());
    }

}
