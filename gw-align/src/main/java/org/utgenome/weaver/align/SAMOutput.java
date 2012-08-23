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
// SAMOutput.java
// Since: 2011/08/16
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.io.OutputStream;
import java.io.PrintWriter;

import org.utgenome.weaver.align.record.AlignmentRecord;
import org.utgenome.weaver.parallel.Reporter;
import org.xerial.lens.SilkLens;
import org.xerial.util.ObjectHandler;
import org.xerial.util.log.Logger;

/**
 * {@link AlignmentRecord} to SAM format converter
 * 
 * @author leo
 * 
 */
public class SAMOutput implements ObjectHandler<AlignmentRecord>, Reporter
{
    private static Logger _logger = Logger.getLogger(SAMOutput.class);

    FMIndexOnGenome       fmIndex;
    PrintWriter           out;
    SequenceBoundary      boundary;
    int                   count   = 0;

    public SAMOutput(SequenceBoundary sequenceBoundary, OutputStream out) {
        this.boundary = sequenceBoundary;
        this.out = new PrintWriter(out, true);
    }

    @Override
    public void init() throws Exception {
        out.print(boundary.toSAMHeader());
    }

    @Override
    public void handle(AlignmentRecord r) throws Exception {
        out.println(r.toSAMLine());
    }

    @Override
    public void finish() throws Exception {
        out.flush();
        out.close();
    }

    @Override
    public void emit(Object result) throws Exception {

        if (_logger.isTraceEnabled())
            _logger.trace(SilkLens.toSilk("result", result));

        if (result != null && result.getClass().isAssignableFrom(AlignmentRecord.class)) {
            AlignmentRecord r = (AlignmentRecord) result;
            out.println(r.toSAMLine());
        }
    }
}
