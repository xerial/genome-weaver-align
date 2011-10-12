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
// KMPMatcher.java
// Since: 2011/10/12
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.util.ArrayList;
import java.util.List;

public class KMPMatcher
{
    public final CharSequence pattern;
    private final int         m;
    private final int[]       prefixIndex;

    public KMPMatcher(CharSequence pattern) {
        this.pattern = pattern;
        this.m = pattern.length();
        prefixIndex = new int[m];
        computePrefixFunction();
    }

    public static List<Integer> find(CharSequence text, CharSequence pattern) {
        return new KMPMatcher(pattern).match(text);
    }

    public List<Integer> match(CharSequence text) {
        final int n = text.length();
        int q = 0;

        List<Integer> matchPos = new ArrayList<Integer>();

        for (int i = 0; i < n; ++i) {
            while (q > 0 && pattern.charAt(q) != text.charAt(i)) {
                q = prefixIndex[q - 1];
            }
            if (pattern.charAt(q) == text.charAt(i))
                ++q;

            if (q >= m) {
                matchPos.add(i - m + 1);
                q = prefixIndex[q - 1];
            }

        }
        return matchPos;
    }

    private void computePrefixFunction() {
        prefixIndex[0] = 0;
        int k = 0;
        for (int i = 1; i < m; ++i) {
            while (k > 0 && pattern.charAt(k) != pattern.charAt(i)) {
                k = prefixIndex[k - 1];
            }

            if (pattern.charAt(k) == pattern.charAt(i))
                ++k;
            prefixIndex[i] = k;
        }
    }

}
