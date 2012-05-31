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
// VariationData.java
// Since: 2011/09/07
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.assembly;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class VariationData
{
    public int                   id;
    public String                chr;
    public int                   start;
    public int                   end;
    public String                refAllele;
    public String                rs;                                     // dbSNP
    public String                nonRefAllele;
    public int                   count;
    public List<AlleleFrequency> freq = new ArrayList<AlleleFrequency>();

    public boolean isSinglePointMutation() {
        return refAllele.length() == 1 && nonRefAllele.length() == 1;
    }

    public static class AlleleFrequency
    {
        public String allele;
        public int    count;
    }

    public static List<VariationData> parse(File file) throws IOException {
        return parse(new BufferedReader(new FileReader(file)));

    }

    public static List<VariationData> parse(BufferedReader r) throws IOException {
        List<VariationData> l = new ArrayList<VariationData>();
        try {
            for (String line; (line = r.readLine()) != null;) {
                String[] s = line.split("\t");
                if (s.length < 7)
                    continue;
                VariationData v = new VariationData();
                v.id = Integer.parseInt(s[0]);
                v.chr = s[1];
                v.start = Integer.parseInt(s[2]);
                v.end = Integer.parseInt(s[3]);
                v.refAllele = s[4];
                v.rs = s[5];
                v.nonRefAllele = s[6];
                for (int i = 7; i < s.length; ++i) {
                    String f = s[i];
                    AlleleFrequency freq = new AlleleFrequency();
                    int cursor = 0;
                    if (f.charAt(cursor++) != '(')
                        break;
                    int nextCursor = cursor + 1;
                    while (nextCursor < f.length() && f.charAt(nextCursor++) != ',') {}
                    freq.allele = f.substring(cursor, nextCursor - 1);

                    cursor = nextCursor;
                    while (nextCursor < f.length() && f.charAt(nextCursor++) != ')') {}
                    freq.count = Integer.parseInt(f.substring(cursor, nextCursor - 1));
                    v.freq.add(freq);
                }
                l.add(v);
            }
        }
        finally {
            r.close();
        }
        return l;
    }

}
