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
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

import org.utgenome.UTGBException;
import org.utgenome.format.fasta.FASTAPullParser;
import org.utgenome.format.fasta.FASTASequence;
import org.utgenome.format.fastq.FastqRead;
import org.utgenome.format.fastq.FastqReader;
import org.xerial.snappy.SnappyInputStream;
import org.xerial.util.ObjectHandler;

public class ReadSequenceReaderFactory
{
    public static ReadSequenceReader createReader(String inputFile) throws IOException, UTGBException {

        boolean gzipped = inputFile.endsWith(".gz");
        boolean snapped = inputFile.endsWith(".snap");

        String prefix = inputFile;
        InputStream in = new FileInputStream(inputFile);
        if (gzipped) {
            prefix = inputFile.replace(".gz$", "");
            in = new GZIPInputStream(in);
        }
        if (snapped) {
            prefix = inputFile.replace(".snap$", "");
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

    public static ReadSequenceReader createFASTAReader(BufferedReader input) {
        return new FASTAReadReader(input);
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

    }

}
