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
// LStringSeq.java
// Since: 2011/02/24
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.sais;

import org.utgenome.weaver.align.LSeq;

public class LStringSeq implements LSeq
{
    private final String s;

    public LStringSeq(String s) {
        this.s = s;
    }

    @Override
    public long lookup(long index) {
        int pos = (int) index;
        if (pos == s.length())
            return 0;
        return s.charAt((int) index);
    }

    @Override
    public long textSize() {
        return s.length() + 1; // add 1 for $ (sentinel)
    }

    @Override
    public void set(long index, long value) {
        throw new UnsupportedOperationException("set");
    }

    @Override
    public long update(long index, long value) {
        throw new UnsupportedOperationException("update");
    }

    @Override
    public String toString() {
        return s + "$";
    }

}
