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

import org.xerial.util.StringUtil;

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
            return String.format("[%s]", StringUtil.join(siF, ","));
        }
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
            return String.format("[%s]", StringUtil.join(siB, ","));
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
            return String.format("[%s]:[%s]", StringUtil.join(siF, ","), StringUtil.join(siB, ","));
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

    }

    public abstract SuffixInterval getForward(ACGT ch);

    public abstract SuffixInterval getBackward(ACGT ch);

    public abstract boolean isEmpty(ACGT ch);

}
