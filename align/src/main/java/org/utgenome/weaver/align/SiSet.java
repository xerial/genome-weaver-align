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
// SiSet.java
// Since: 2011/08/12
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

/**
 * Various types of suffix intervals for forward, backward, bidirectional
 * searches.
 * 
 * @author leo
 * 
 */
public abstract class SiSet
{
    public static class ForwardSiSet extends SiSet
    {
        private final SuffixInterval[] siF;

        public ForwardSiSet(SuffixInterval[] siF) {
            this.siF = siF;
        }

        @Override
        public boolean isEmpty(ACGT ch) {
            return (siF[ch.code] == null) || siF[ch.code].isEmpty();
        }

        @Override
        public SuffixInterval getForward(ACGT ch) {
            return siF[ch.code];
        }

        @Override
        public SuffixInterval getBackward(ACGT ch) {
            return null;
        }

        @Override
        public String toString() {
            return "F" + toString(siF);
        }

        @Override
        public SuffixInterval getNext(ACGT ch) {
            return getForward(ch);
        }
    }

    protected static String toString(SuffixInterval[] si) {
        StringBuilder s = new StringBuilder();
        s.append("[");
        int count = 0;
        for (int i = 0; i < si.length; ++i) {
            if (si[i] != null) {
                if (count != 0)
                    s.append(" ");
                s.append(ACGT.decode(i));
                s.append(":");
                s.append(si[i].toString());
                ++count;
            }
        }
        s.append("]");
        return s.toString();
    }

    public static class BackwardSiSet extends SiSet
    {
        private final SuffixInterval[] siB;

        public BackwardSiSet(SuffixInterval[] siB) {
            this.siB = siB;
        }

        @Override
        public SuffixInterval getForward(ACGT ch) {
            return null;
        }

        @Override
        public SuffixInterval getBackward(ACGT ch) {
            return siB[ch.code];
        }

        @Override
        public boolean isEmpty(ACGT ch) {
            return (siB[ch.code] == null) || siB[ch.code].isEmpty();
        }

        @Override
        public String toString() {
            return "B" + toString(siB);
        }

        @Override
        public SuffixInterval getNext(ACGT ch) {
            return getBackward(ch);
        }
    }

    public static class BidirectionalSiSet extends SiSet
    {
        private final SuffixInterval[] siF;
        private final SuffixInterval[] siB;

        public BidirectionalSiSet(SuffixInterval[] siF, SuffixInterval[] siB) {
            this.siF = siF;
            this.siB = siB;
        }

        @Override
        public SuffixInterval getForward(ACGT ch) {
            return siF[ch.code];
        }

        @Override
        public SuffixInterval getBackward(ACGT ch) {
            return siB[ch.code];
        }

        @Override
        public boolean isEmpty(ACGT ch) {
            return (siF[ch.code] == null) || siF[ch.code].isEmpty();
        }

        @Override
        public String toString() {
            return String.format("F%s:B%s", toString(siF), toString(siB));
        }

        @Override
        public SuffixInterval getNext(ACGT ch) {
            return getForward(ch);
        }
    }

    public static SiSet empty = new EmptySiSet();

    public static class EmptySiSet extends SiSet
    {

        @Override
        public SuffixInterval getForward(ACGT ch) {
            return null;
        }

        @Override
        public SuffixInterval getBackward(ACGT ch) {
            return null;
        }

        @Override
        public boolean isEmpty(ACGT ch) {
            return true;
        }

        @Override
        public String toString() {
            return "[]";
        }

        @Override
        public SuffixInterval getNext(ACGT ch) {
            return null;
        }

    }

    public abstract SuffixInterval getNext(ACGT ch);

    public abstract SuffixInterval getForward(ACGT ch);

    public abstract SuffixInterval getBackward(ACGT ch);

    public abstract boolean isEmpty(ACGT ch);

}
