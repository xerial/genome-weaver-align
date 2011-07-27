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
// IUPACBinaryInfo.java
// Since: 2011/02/14
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

import org.utgenome.UTGBErrorCode;
import org.utgenome.UTGBException;
import org.xerial.lens.SilkLens;

/**
 * Holding sequence boundaries of concatenated sequences
 * 
 * @author leo
 * 
 */
public class SequenceBoundary
{
    public static class SequenceIndex
    {
        public String name;
        public String desc;
        public long   length;
        public long   offset;

        public SequenceIndex(String name, String desc, long length, long offset) {
            this.name = name;
            this.desc = desc;
            this.length = length;
            this.offset = offset;
        }
    }

    public List<SequenceIndex>    index;
    public long                   totalSize;

    private TreeMap<Long, String> indexToChrTable;

    public static class PosOnGenome
    {
        public final String chr;
        public final int    pos;
        public final Strand strand;

        public PosOnGenome(String chr, int pos, Strand strand) {
            this.chr = chr;
            this.pos = pos;
            this.strand = strand;
        }
    }

    public String toSAMHeader() {
        StringWriter buf = new StringWriter();
        for (SequenceIndex each : index) {
            buf.append(String.format("@SQ\tSN:%s\tLN:%d\n", each.name, each.length));
        }
        return buf.toString();
    }

    /**
     * @param textIndex
     * @return 1-origin index
     * @throws UTGBException
     */
    public PosOnGenome translate(long textIndex, Strand strand) throws UTGBException {

        if (indexToChrTable == null) {
            indexToChrTable = new TreeMap<Long, String>();
            for (SequenceIndex each : index) {
                indexToChrTable.put(each.offset, each.name);
            }
        }

        SortedMap<Long, String> headMap = indexToChrTable.headMap(textIndex);
        if (headMap == null || headMap.isEmpty())
            throw new UTGBException(UTGBErrorCode.INVALID_INPUT, "invalid index: " + textIndex);

        long offset = headMap.lastKey();
        String chr = headMap.get(offset);
        int start = (int) (textIndex - offset);
        return new PosOnGenome(chr, start, strand);
    }

    public static SequenceBoundary loadSilk(File silkFile) throws UTGBException {
        BufferedReader input = null;
        try {
            try {
                input = new BufferedReader(new FileReader(silkFile));
                return SilkLens.loadSilk(SequenceBoundary.class, input);
            }
            finally {
                if (input != null)
                    input.close();
            }
        }
        catch (Exception e) {
            throw UTGBException.convert(e);
        }

    }

    public static SequenceBoundary createFromSingleSeq(String name, ACGTSequence seq) {
        SequenceBoundary s = new SequenceBoundary();
        s.totalSize = seq.textSize();
        s.index = new ArrayList<SequenceBoundary.SequenceIndex>();
        s.index.add(new SequenceIndex(name, name, seq.textSize(), 0));
        return s;
    }

}
