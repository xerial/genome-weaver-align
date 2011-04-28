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

import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;

import org.utgenome.UTGBException;
import org.utgenome.util.StandardOutputStream;
import org.utgenome.weaver.align.record.AlignmentRecord;
import org.utgenome.weaver.align.record.RawRead;
import org.utgenome.weaver.align.record.ReadSequenceReader;
import org.utgenome.weaver.align.record.ReadSequenceReaderFactory;
import org.utgenome.weaver.align.strategy.BWAStrategy;
import org.xerial.util.ObjectHandler;
import org.xerial.util.ObjectHandlerBase;
import org.xerial.util.log.Logger;
import org.xerial.util.opt.Argument;
import org.xerial.util.opt.Option;

/**
 * Burrows-Wheeler aligner for IUPAC sequences
 * 
 * @author leo
 * 
 */
public class BWAlign extends GenomeWeaverCommand
{
    static Logger _logger = Logger.getLogger(BWAlign.class);

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

    @Argument(index = 1)
    private String readFile;

    @Option(symbol = "q", description = "query sequence")
    private String query;

    //    @Option(longName = "sam", description = "output in SAM format")
    //    public boolean outputSAM           = false;

    @Option(symbol = "N", description = "Num mismatches allowed. default=0")
    public int     numMismachesAllowed = 0;

    public static class SAMOutput implements ObjectHandler<AlignmentRecord>
    {

        FMIndexOnGenome fmIndex;
        PrintWriter     out;

        public SAMOutput(FMIndexOnGenome fmIndex, OutputStream out) {
            this.fmIndex = fmIndex;
            this.out = new PrintWriter(new OutputStreamWriter(out));
        }

        @Override
        public void init() throws Exception {
            fmIndex.outputSAMHeader(out);
        }

        @Override
        public void handle(AlignmentRecord r) throws Exception {
            out.println(r.toSAMLine());
        }

        @Override
        public void finish() throws Exception {
            out.close();
        }
    }

    @Override
    public void execute(String[] args) throws Exception {

        if (query == null && readFile == null)
            throw new UTGBException("no query is given");

        ReadSequenceReader reader = null;
        if (query != null) {
            _logger.info("query sequence: " + query);
            reader = ReadSequenceReaderFactory.singleQueryReader(query);
        }
        else if (readFile != null) {
            reader = ReadSequenceReaderFactory.createReader(readFile);
        }

        FMIndexOnGenome fmIndex = new FMIndexOnGenome(fastaFilePrefix);
        SAMOutput samOutput = new SAMOutput(fmIndex, new StandardOutputStream());
        query(fmIndex, reader, samOutput);
    }

    private static class GenomeCoordinateConverter extends ObjectHandlerBase<AlignmentSA>
    {

        private FMIndexOnGenome                fmIndex;
        private ObjectHandler<AlignmentRecord> handler;

        public GenomeCoordinateConverter(FMIndexOnGenome fmIndex, ObjectHandler<AlignmentRecord> handler) {
            this.fmIndex = fmIndex;
            this.handler = handler;
        }

        @Override
        public void handle(AlignmentSA aln) throws Exception {
            fmIndex.toGenomeCoordinate(aln, handler);
        }

    }

    public static void query(final FMIndexOnGenome fmIndex, ReadSequenceReader readReader,
            final ObjectHandler<AlignmentRecord> handler) throws Exception {

        handler.init();
        final BWAStrategy aligner = new BWAStrategy(fmIndex);
        readReader.parse(new ObjectHandlerBase<RawRead>() {
            @Override
            public void handle(final RawRead input) throws Exception {
                aligner.align(input, new GenomeCoordinateConverter(fmIndex, handler));
            };
        });
        handler.finish();

    }

    public static void querySingle(String fastaFilePrefix, final String query,
            final ObjectHandler<AlignmentRecord> resultHandler) throws Exception {

        final FMIndexOnGenome fmIndex = new FMIndexOnGenome(fastaFilePrefix);
        query(fmIndex, ReadSequenceReaderFactory.singleQueryReader(query), resultHandler);
    }

}
