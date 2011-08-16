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
// ReadSequenceReaderFactory.java
// Since: 2011/04/28
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.record;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

import org.utgenome.UTGBErrorCode;
import org.utgenome.UTGBException;
import org.utgenome.format.fasta.FASTAPullParser;
import org.utgenome.format.fasta.FASTASequence;
import org.utgenome.format.fastq.FastqRead;
import org.utgenome.format.fastq.FastqReader;
import org.xerial.snappy.SnappyInputStream;
import org.xerial.util.ObjectHandler;

public class ReadSequenceReaderFactory
{
    public static ReadSequenceReader createReader(String[] inputFiles) throws IOException, UTGBException {
        switch (inputFiles.length) {
        case 1:
            return createReader(inputFiles[0]);
        case 2:
            return createPEReader(inputFiles);
        default:
            throw new UTGBException(UTGBErrorCode.INVALID_INPUT,
                    "# of input read files must be one (single-end) or two (paired-end).");
        }
    }

    private static ReadSequenceReader createPEReader(String[] inputFiles) throws IOException, UTGBException {

        assert (inputFiles.length == 2);

        ReadSequenceReader r1 = createReader(inputFiles[0]);
        ReadSequenceReader r2 = createReader(inputFiles[1]);

        return new ReadSequenceReader() {
            @Override
            public RawRead next() throws Exception {

                // TODO Auto-generated method stub
                return null;
            }

            @Override
            public void parse(ObjectHandler<RawRead> handler) throws Exception {
                // TODO Auto-generated method stub

            }

            @Override
            public void close() throws IOException {
                // TODO Auto-generated method stub

            }
        };
    }

    public static class PairedReadReader implements ReadSequenceReader
    {
        private final ReadSequenceReader p1;
        private final ReadSequenceReader p2;

        public PairedReadReader(ReadSequenceReader p1, ReadSequenceReader p2) {
            this.p1 = p1;
            this.p2 = p2;
        }

        @Override
        public RawRead next() throws Exception {
            RawRead r1 = p1.next();
            RawRead r2 = p2.next();
            if (r1 != null && r2 != null)
                return new PairedEndRead(r1, r2);
            else
                return null;
        }

        @Override
        public void parse(ObjectHandler<RawRead> handler) throws Exception {
            handler.init();
            for (RawRead r1, r2; (r1 = p1.next()) != null && (r2 = p2.next()) != null;) {
                handler.handle(new PairedEndRead(r1, r2));
            }
            handler.finish();
        }

        @Override
        public void close() throws IOException {
            p1.close();
            p2.close();
        }

    }

    public static ReadSequenceReader createReader(String inputFile) throws IOException, UTGBException {

        boolean gzipped = inputFile.endsWith(".gz");
        boolean snapped = inputFile.endsWith(".snap");

        String prefix = inputFile;
        InputStream in = new FileInputStream(inputFile);
        if (gzipped) {
            prefix = inputFile.replaceAll("\\.gz$", "");
            in = new GZIPInputStream(in);
        }
        if (snapped) {
            prefix = inputFile.replaceAll("\\.snap$", "");
            in = new SnappyInputStream(in);
        }

        BufferedReader reader = new BufferedReader(new InputStreamReader(in), 4 * 1024 * 1024); // Use 4MB buffer

        if (prefix.endsWith(".fa") || prefix.endsWith(".fasta") || prefix.endsWith(".fan")) {
            // FASTA file
            return new FASTAReadReader(reader);
        }
        else if (prefix.endsWith(".fastq") || prefix.endsWith(".fq")) {
            return createFASTQReader(reader);
        }

        throw new UTGBException("Unsupported file type: " + inputFile);
    }

    public static ReadSequenceReader createFASTQReader(BufferedReader input) {
        return new FASTQReadReader(input);
    }

    public static ReadSequenceReader createFASTQReader(String fastqFile) throws FileNotFoundException {
        return new FASTQReadReader(new BufferedReader(new FileReader(fastqFile)));
    }

    public static ReadSequenceReader createFASTAReader(BufferedReader input) {
        return new FASTAReadReader(input);
    }

    public static ReadSequenceReader singleQueryReader(final String query) {
        return new ReadSequenceReader() {
            int count = 0;

            @Override
            public void parse(ObjectHandler<RawRead> handler) throws Exception {
                handler.handle(new ReadSequence(query, query, null));
                count++;
            }

            @Override
            public void close() throws IOException {
                // do nothing
            }

            @Override
            public RawRead next() throws Exception {
                if (count++ == 0)
                    return new ReadSequence(query, query, null);
                return null;
            }
        };
    }

    private static class FASTAReadReader implements ReadSequenceReader
    {
        FASTAPullParser input;

        public FASTAReadReader(BufferedReader input) {
            this.input = new FASTAPullParser(input);
        }

        @Override
        public void parse(ObjectHandler<RawRead> handler) throws Exception {
            handler.init();
            for (FASTASequence seq; (seq = input.nextSequence()) != null;) {
                ReadSequence read = ReadSequence.createFrom(seq);
                handler.handle(read);
            }
            handler.finish();
        }

        @Override
        public void close() throws IOException {
            input.close();
        }

        @Override
        public RawRead next() throws Exception {
            FASTASequence seq = input.nextSequence();
            if (seq == null)
                return null;
            else
                return ReadSequence.createFrom(seq);
        }
    }

    private static class FASTQReadReader implements ReadSequenceReader
    {
        FastqReader reader;

        public FASTQReadReader(BufferedReader in) {
            reader = new FastqReader(in);
        }

        @Override
        public void parse(ObjectHandler<RawRead> handler) throws Exception {
            handler.init();
            for (FastqRead read; (read = reader.next()) != null;) {
                handler.handle(ReadSequence.createFrom(read));
            }
            handler.finish();
        }

        @Override
        public void close() throws IOException {
            reader.close();
        }

        @Override
        public RawRead next() throws Exception {
            FastqRead read = reader.next();
            if (read != null)
                return ReadSequence.createFrom(read);
            else
                return null;
        }

    }

}
