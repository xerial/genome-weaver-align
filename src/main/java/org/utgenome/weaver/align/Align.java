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
import org.utgenome.weaver.align.strategy.BidirectionalSuffixFilter;
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
        return "read alignment";
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

        OutputStream out = config.silent ? new NullOutputStream() : new StandardOutputStream();
        SAMOutput reporter = new SAMOutput(SequenceBoundary.load(config.refSeq), out);
        CommonDataSet common = CommonDataSet.prepare(config, reporter);

        try {
            reporter.init();
            query(common, reader);
        }
        finally {
            reporter.finish();
        }
    }

    public static void querySingle(FMIndexOnGenome fmIndex, ACGTSequence reference, String query, Reporter out)
            throws Exception {
        query(new CommonDataSet(fmIndex, reference, new AlignmentConfig(), out),
                ReadReaderFactory.singleQueryReader(query));
    }

    public static void query(CommonDataSet common, ReadReader readReader) throws Exception {

        _logger.debug("Alignment mode: %s", common.config.strategy.description);
        Aligner aligner = null;
        switch (common.config.strategy) {
        case BSF:
            aligner = new BidirectionalSuffixFilter(common.fmIndex, common.reference, common.config);
            break;
        case SF:
            aligner = new SuffixFilter(common.fmIndex, common.reference, common.config);
            break;
        case BD:
            aligner = new BidirectionalBWT(common.fmIndex, common.reporter);
            break;
        case BWA:
            aligner = new BWAAligner(fmIndex, config, reporter);
            ((BWAAligner) aligner).aligner.disableBidirectionalSearch();
            break;
        default:
            throw new UTGBException(String.format("%s mode is not supported", common.config.strategy));
            //        case BD:
            //            aligner = new BWAAligner(fmIndex, config, reporter);
            //            break;
        }

        readReader.parse(new PassReadToAligner(common, aligner));
    }

    public static class CommonDataSet
    {
        private final FMIndexOnGenome fmIndex;
        private final ACGTSequence    reference;
        private final AlignmentConfig config;
        private Reporter              reporter;

        private int                   count = 0;
        private StopWatch             timer = new StopWatch();

        public CommonDataSet(FMIndexOnGenome fmIndex, ACGTSequence reference, AlignmentConfig config, Reporter reporter) {
            this.fmIndex = fmIndex;
            this.reference = reference;
            this.config = config;
            this.reporter = reporter;
        }

        public static CommonDataSet prepare(AlignmentConfig config, Reporter out) throws Exception {
            BWTFiles forwardDB = new BWTFiles(config.refSeq, Strand.FORWARD);
            SequenceBoundary b = SequenceBoundary.loadSilk(forwardDB.pacIndex());

            FMIndexOnGenome fmIndex = FMIndexOnGenome.load(config.refSeq);

            _logger.info("loading reference sequence %s", forwardDB.pac());
            ACGTSequence reference = ACGTSequence.loadFrom(forwardDB.pac());

            return new CommonDataSet(fmIndex, reference, config, out);

        }

    }

    protected static class PassReadToAligner extends ObjectHandlerBase<Read>
    {
        private final CommonDataSet common;
        private final Aligner       aligner;

        private int                 count = 0;
        private StopWatch           timer = new StopWatch();

        public PassReadToAligner(CommonDataSet common, Aligner aligner) {
            this.common = common;
            this.aligner = aligner;
        }

        @Override
        public void handle(Read input) throws Exception {
            aligner.align(input, common.reporter);
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

}
