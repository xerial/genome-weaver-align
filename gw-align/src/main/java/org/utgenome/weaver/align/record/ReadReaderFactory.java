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
import java.util.List;
import java.util.zip.GZIPInputStream;

import org.utgenome.UTGBErrorCode;
import org.utgenome.UTGBException;
import org.utgenome.format.fasta.FASTAPullParser;
import org.utgenome.format.fasta.FASTASequence;
import org.utgenome.format.fastq.FastqRead;
import org.utgenome.format.fastq.FastqReader;
import org.xerial.snappy.SnappyInputStream;
import org.xerial.util.ObjectHandler;

public class ReadReaderFactory
{
    public static ReadReader createReader(List<String> inputFiles) throws IOException, UTGBException {
        switch (inputFiles.size()) {
        case 1:
            return createReader(inputFiles.get(0));
        case 2:
            return createPEReader(inputFiles);
        default:
            throw new UTGBException(UTGBErrorCode.INVALID_INPUT,
                    "# of input read files must be one (single-end) or two (paired-end).");
        }
    }

    private static ReadReader createPEReader(List<String> inputFiles) throws IOException, UTGBException {

        assert (inputFiles.size() == 2);

        ReadReader r1 = createReader(inputFiles.get(0));
        ReadReader r2 = createReader(inputFiles.get(1));

        return new ReadReader() {
            @Override
            public Read next() throws Exception {

                // TODO Auto-generated method stub
                return null;
            }

            @Override
            public void parse(ObjectHandler<Read> handler) throws Exception {
                // TODO Auto-generated method stub

            }

            @Override
            public void close() throws IOException {
                // TODO Auto-generated method stub

            }
        };
    }

    public static class PairedReadReader implements ReadReader
    {
        private final ReadReader p1;
        private final ReadReader p2;

        public PairedReadReader(ReadReader p1, ReadReader p2) {
            this.p1 = p1;
            this.p2 = p2;
        }

        @Override
        public Read next() throws Exception {
            Read r1 = p1.next();
            Read r2 = p2.next();
            if (r1 != null && r2 != null)
                return new PairedEndRead(r1, r2);
            else
                return null;
        }

        @Override
        public void parse(ObjectHandler<Read> handler) throws Exception {
            handler.init();
            for (Read r1, r2; (r1 = p1.next()) != null && (r2 = p2.next()) != null;) {
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

    public static ReadReader createReader(String inputFile) throws IOException, UTGBException {

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

    public static ReadReader createFASTQReader(BufferedReader input) {
        return new FASTQReadReader(input);
    }

    public static ReadReader createFASTQReader(String fastqFile) throws FileNotFoundException {
        return new FASTQReadReader(new BufferedReader(new FileReader(fastqFile)));
    }

    public static ReadReader createFASTAReader(BufferedReader input) {
        return new FASTAReadReader(input);
    }

    public static ReadReader singleQueryReader(final String query) {
        return new ReadReader() {
            int count = 0;

            @Override
            public void parse(ObjectHandler<Read> handler) throws Exception {
                handler.init();
                handler.handle(new SingleEndRead("read", query, null));
                count++;
                handler.finish();
            }

            @Override
            public void close() throws IOException {
                // do nothing
            }

            @Override
            public Read next() throws Exception {
                if (count++ == 0)
                    return new SingleEndRead(query, query, null);
                return null;
            }
        };
    }

    private static class FASTAReadReader implements ReadReader
    {
        FASTAPullParser input;

        public FASTAReadReader(BufferedReader input) {
            this.input = new FASTAPullParser(input);
        }

        @Override
        public void parse(ObjectHandler<Read> handler) throws Exception {
            handler.init();
            for (FASTASequence seq; (seq = input.nextSequence()) != null;) {
                SingleEndRead read = SingleEndRead.createFrom(seq);
                handler.handle(read);
            }
            handler.finish();
        }

        @Override
        public void close() throws IOException {
            input.close();
        }

        @Override
        public Read next() throws Exception {
            FASTASequence seq = input.nextSequence();
            if (seq == null)
                return null;
            else
                return SingleEndRead.createFrom(seq);
        }
    }

    private static class FASTQReadReader implements ReadReader
    {
        FastqReader reader;

        public FASTQReadReader(BufferedReader in) {
            reader = new FastqReader(in);
        }

        @Override
        public void parse(ObjectHandler<Read> handler) throws Exception {
            handler.init();
            for (FastqRead read; (read = reader.next()) != null;) {
                handler.handle(SingleEndRead.createFrom(read));
            }
            handler.finish();
        }

        @Override
        public void close() throws IOException {
            reader.close();
        }

        @Override
        public Read next() throws Exception {
            FastqRead read = reader.next();
            if (read != null)
                return SingleEndRead.createFrom(read);
            else
                return null;
        }

    }

}
