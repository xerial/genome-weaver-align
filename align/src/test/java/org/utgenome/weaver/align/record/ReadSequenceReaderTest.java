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
// ReadSequenceReaderTest.java
// Since: 2011/04/28
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.record;

import static org.junit.Assert.*;

import org.junit.Test;
import org.xerial.util.FileResource;
import org.xerial.util.ObjectHandlerBase;
import org.xerial.util.log.Logger;

public class ReadSequenceReaderTest
{
    private static Logger _logger = Logger.getLogger(ReadSequenceReaderTest.class);

    @Test
    public void read() throws Exception {
        ReadReader r = ReadReaderFactory.createFASTQReader(FileResource.open(
                ReadSequenceReaderTest.class, "sample.fastq"));
        r.parse(new ObjectHandlerBase<Read>() {
            int readCount = 0;

            @Override
            public void handle(Read input) throws Exception {
                _logger.info(input);
                readCount++;
            }

            @Override
            public void finish() throws Exception {
                assertEquals(3, readCount);
            }
        });
    }

    @Test
    public void readFASTA() throws Exception {
        ReadReader r = ReadReaderFactory.createFASTAReader(FileResource.open(
                ReadSequenceReaderTest.class, "sample.fa"));
        r.parse(new ObjectHandlerBase<Read>() {
            int readCount = 0;

            @Override
            public void handle(Read input) throws Exception {
                _logger.info(input);
                readCount++;
            }

            @Override
            public void finish() throws Exception {
                assertEquals(2, readCount);
            }
        });
    }
}
