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
// StringSeq.java
// Since: 2011/04/07
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

/**
 * String wrapper
 * 
 * @author leo
 * 
 */
public class StringSeq implements LSeq
{
    private String s;

    public StringSeq(String s) {
        this.s = s;
    }

    @Override
    public String toString() {
        return s;
    }

    @Override
    public long lookup(long index) {
        return s.charAt((int) index);
    }

    @Override
    public long textSize() {
        return s.length();
    }

    @Override
    public void set(long index, long value) {
        throw new UnsupportedOperationException("set");
    }

    @Override
    public long increment(long index, long value) {
        throw new UnsupportedOperationException("update");
    }

}
