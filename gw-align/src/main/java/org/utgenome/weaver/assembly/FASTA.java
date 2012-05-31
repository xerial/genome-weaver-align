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
// FASTA.java
// Since: 2011/10/05
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.assembly;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;

import org.utgenome.format.fasta.CompactFASTA;
import org.utgenome.format.fasta.FASTAPullParser;
import org.utgenome.weaver.align.ACGTSequence;

/**
 * FASTA sequence file in memory
 * 
 * @author leo
 * 
 */
public class FASTA
{
    private final LinkedHashMap<String, ACGTSequence> sequenceTable;

    private FASTA(LinkedHashMap<String, ACGTSequence> sequenceTable) {
        this.sequenceTable = sequenceTable;
    }

    public ACGTSequence getACGTSequence(String chr) {
        return sequenceTable.get(chr);
    }

    public void toFASTA(Writer out) throws IOException {
        List<String> chrList = new ArrayList<String>(sequenceTable.keySet());
        Collections.sort(chrList, new Comparator4ChrName());

        final int width = 50;

        for (String chr : chrList) {
            ACGTSequence seq = sequenceTable.get(chr);
            out.append(String.format(">%s\n", chr));
            for (long i = 0; i < seq.length(); i += width) {
                long seqEnd = Math.min(i + width, seq.length());
                ACGTSequence line = seq.subString(i, seqEnd);
                out.append(line.toString());
                out.append("\n");
            }
        }
        out.flush();
    }

    public static FASTA load(File fastaFile) throws IOException {
        LinkedHashMap<String, ACGTSequence> table = new LinkedHashMap<String, ACGTSequence>();
        FASTAPullParser fasta = new FASTAPullParser(fastaFile);
        int lineCount = 0;
        for (String desc; (desc = fasta.nextDescriptionLine()) != null; lineCount++) {
            String chr = CompactFASTA.pickSequenceName(desc);
            ACGTSequence seq = new ACGTSequence();
            for (String seqLine; (seqLine = fasta.nextSequenceLine()) != null; lineCount++) {
                seqLine = seqLine.trim();
                for (int i = 0; i < seqLine.length(); ++i) {
                    seq.append(seqLine.charAt(i));
                }
            }
            table.put(chr, seq);
        }

        return new FASTA(table);
    }

    public static class Comparator4ChrName implements Comparator<String>
    {

        public static int compareChrName(String p, String q) {
            int x = p.length();
            int y = q.length();
            int n = Math.min(x, y);
            for (int i = 0; i < n; i++) {
                char c = p.charAt(i);
                char d = q.charAt(i);
                if (c != d) {
                    boolean f = (c >= '0' && c <= '9');
                    boolean g = (d >= '0' && d <= '9');
                    if (f && !g) {
                        return -1;
                    }
                    else if (!f && g) {
                        return 1;
                    }
                    else if (!f && !g) {
                        return c - d;
                    }
                    if (x != y) {
                        return x - y;
                    }
                    return c - d;
                }
            }
            return x - y;
        }

        public int compare(String a, String b) {
            return compareChrName(a, b);
        }
    }

}
