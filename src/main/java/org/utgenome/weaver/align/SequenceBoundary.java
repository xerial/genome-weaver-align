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
    public int                    totalSize;

    private TreeMap<Long, String> indexToChrTable;

    public static File getFileName(String prefix) {
        return new File(prefix + ".i.silk");
    }

    public static class PosOnGenome
    {
        public final String chr;
        public final int    pos;

        public PosOnGenome(String chr, int pos) {
            this.chr = chr;
            this.pos = pos;
        }
    }

    /**
     * @param textIndex
     * @return 1-origin index
     * @throws UTGBException
     */
    public PosOnGenome translate(long textIndex, boolean isReverse) throws UTGBException {

        if (indexToChrTable == null) {
            indexToChrTable = new TreeMap<Long, String>();
            for (SequenceIndex each : index) {
                indexToChrTable.put(each.offset, each.name);
            }
        }

        if (isReverse)
            textIndex = totalSize - textIndex - 1;

        SortedMap<Long, String> headMap = indexToChrTable.headMap(textIndex + 1);
        if (headMap == null || headMap.isEmpty())
            throw new UTGBException(UTGBErrorCode.INVALID_INPUT, "invalid index: " + textIndex);

        long offset = headMap.lastKey();
        String chr = headMap.get(offset);
        int start = (int) (textIndex - offset + 1);
        return new PosOnGenome(chr, start);
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

}
