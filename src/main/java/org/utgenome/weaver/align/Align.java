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
// SFAlign.java
// Since: 2011/08/05
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.io.OutputStream;

import org.utgenome.UTGBException;
import org.utgenome.weaver.GenomeWeaverCommand;
import org.utgenome.weaver.align.record.AlignmentRecord;
import org.utgenome.weaver.align.record.Read;
import org.utgenome.weaver.align.record.ReadReader;
import org.utgenome.weaver.align.record.ReadReaderFactory;
import org.utgenome.weaver.align.strategy.BWAState;
import org.utgenome.weaver.align.strategy.BidirectionalBWT;
import org.utgenome.weaver.align.strategy.SuffixFilter;
import org.utgenome.weaver.parallel.Reporter;
import org.xerial.lens.SilkLens;
import org.xerial.util.ObjectHandler;
import org.xerial.util.ObjectHandlerBase;
import org.xerial.util.StopWatch;
import org.xerial.util.io.NullOutputStream;
import org.xerial.util.io.StandardOutputStream;
import org.xerial.util.log.Logger;

/**
 * Alignment command
 * 
 * @author leo
 * 
 */
public class Align extends GenomeWeaverCommand
{
    private static Logger _logger = Logger.getLogger(Align.class);

    @Override
    public String name() {
        return "align";
    }

    @Override
    public String getOneLineDescription() {
        return "suffix filter alignment";
    }

    @Override
    public Object getOptionHolder() {
        return config;
    }

    private AlignmentConfig config = new AlignmentConfig();

    @Override
    public void execute(String[] args) throws Exception {

        if (config.query == null && config.readFiles == null)
            throw new UTGBException("no query is given");

        ReadReader reader = null;
        if (config.query != null) {
            _logger.info("query sequence: " + config.query);
            reader = ReadReaderFactory.singleQueryReader(config.query);
        }
        else if (config.readFiles != null && !config.readFiles.isEmpty()) {
            reader = ReadReaderFactory.createReader(config.readFiles);
        }
        else {
            throw new UTGBException("no query is given");
        }

        BWTFiles forwardDB = new BWTFiles(config.refSeq, Strand.FORWARD);
        SequenceBoundary b = SequenceBoundary.loadSilk(forwardDB.pacIndex());

        FMIndexOnGenome fmIndex = FMIndexOnGenome.load(config.refSeq);
        OutputStream out = config.quiet ? new NullOutputStream() : new StandardOutputStream();
        SAMOutput reporter = new SAMOutput(fmIndex.getSequenceBoundary(), out);
        try {
            reporter.init();
            _logger.info("loading reference sequence %s", forwardDB.pac());
            ACGTSequence reference = ACGTSequence.loadFrom(forwardDB.pac());
            query(fmIndex, reference, config, reader, reporter);
        }
        finally {
            reporter.finish();
        }
    }

    public static void querySingle(FMIndexOnGenome fmIndex, ACGTSequence reference, String query, Reporter out)
            throws Exception {
        query(fmIndex, reference, new AlignmentConfig(), ReadReaderFactory.singleQueryReader(query), out);
    }

    public static void query(FMIndexOnGenome fmIndex, ACGTSequence reference, AlignmentConfig config,
            ReadReader readReader, Reporter reporter) throws Exception {

        ObjectHandler<Read> aligner = null;
        switch (config.strategy) {
        default:
        case SF:
            aligner = new SuffixFilterAligner(fmIndex, reference, config, reporter);
            break;
        case BWA:
            aligner = new BWAAligner(fmIndex, config, reporter);
            ((BWAAligner) aligner).aligner.disableBidirectionalSearch();
            break;
        case BD:
            aligner = new BWAAligner(fmIndex, config, reporter);
            break;
        }
        _logger.debug("Alignment mode: %s", config.strategy.description);

        readReader.parse(aligner);
    }

    public static class SuffixFilterAligner extends ObjectHandlerBase<Read>
    {

        private final FMIndexOnGenome fmIndex;
        private final ACGTSequence    reference;
        private final AlignmentConfig config;
        private Reporter              reporter;

        private int                   count = 0;
        private StopWatch             timer = new StopWatch();
        private SuffixFilter          sf;

        public SuffixFilterAligner(FMIndexOnGenome fmIndex, ACGTSequence reference, AlignmentConfig config,
                Reporter reporter) {
            this.fmIndex = fmIndex;
            this.reference = reference;
            this.config = config;
            this.reporter = reporter;
        }

        @Override
        public void handle(Read read) throws Exception {
            if (sf == null)
                sf = new SuffixFilter(fmIndex, reference, config);
            sf.align(read, reporter);
            count++;
            double time = timer.getElapsedTime();
            if (count % 10000 == 0) {
                _logger.info("%,d reads are processed in %.2f sec. %,.0f reads/sec.", count, time, count / time);
            }
        }

    }

    private static class GenomeCoordinateConverter extends ObjectHandlerBase<BWAState>
    {
        private FMIndexOnGenome                fmIndex;
        private ObjectHandler<AlignmentRecord> handler;

        public GenomeCoordinateConverter(FMIndexOnGenome fmIndex, ObjectHandler<AlignmentRecord> handler) {
            this.fmIndex = fmIndex;
            this.handler = handler;
        }

        @Override
        public void handle(BWAState aln) throws Exception {
            if (_logger.isTraceEnabled())
                _logger.info(SilkLens.toSilk("alignment", aln));

            aln.toGenomeCoordinate(fmIndex, handler);
        }
    }

    public static class BWAAligner extends ObjectHandlerBase<Read>
    {
        private final FMIndexOnGenome fmIndex;
        private final AlignmentConfig config;
        private Reporter              reporter;

        private int                   count = 0;
        private StopWatch             timer = new StopWatch();

        public final BidirectionalBWT aligner;

        public BWAAligner(FMIndexOnGenome fmIndex, AlignmentConfig config, Reporter reporter) {
            this.fmIndex = fmIndex;
            this.config = config;
            this.reporter = reporter;
            aligner = new BidirectionalBWT(fmIndex, reporter);
        }

        @Override
        public void handle(Read input) throws Exception {
            aligner.align(input);
            count++;
            double time = timer.getElapsedTime();
            if (count % 10000 == 0) {
                _logger.info(String.format("%,d reads are processed in %.2f sec.", count, time));
            }
        }

    }
}
