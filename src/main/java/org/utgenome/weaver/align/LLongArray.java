/*--------------------------------------------------------------------------
 *  Copyright 2010 utgenome.org
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
// utgb-core Project
//
// UInt32Array.java
// Since: 2010/11/05
//
//--------------------------------------
package org.utgenome.weaver.align;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Array capable to store more than 2G (2 x 1024 x 1024 x 1024) entries
 * 
 * @author leo
 * 
 */
public class LLongArray
{

    private final long        size;

    private ArrayList<long[]> array;

    public LLongArray(long size) {
        this.size = size;

        // (flag)|---(array index)----|---------------------------|
        // (flag)|------(31 bit)-------|------(index: 31bit)------|
        int container = (int) (size >> 31);
        int remainder = (int) size & 0x7FFFFFFF;

        array = new ArrayList<long[]>(container + 1);
        for (int i = 0; i < container; i++) {
            array.add(new long[Integer.MAX_VALUE]);
        }
        if (remainder > 0)
            array.add(new long[remainder]);
    }

    public long size() {
        return size;
    }

    private long[] container(long index) {
        return array.get((int) (index >> 31));
    }

    private static int offset(long index) {
        return (int) index & 0x7FFFFFFF;
    }

    public long get(long index) {
        return container(index)[offset(index)];
    }

    public void set(long index, long value) {
        container(index)[offset(index)] = value;
    }

    public void setOR(long index, long value) {
        container(index)[offset(index)] |= value;
    }

    public void setAND(long index, long value) {
        container(index)[offset(index)] &= value;
    }

    public void fill(long value) {
        for (int c = 0; c < array.size(); ++c) {
            long[] container = array.get(c);
            for (int i = 0; i < container.length; ++i) {
                container[i] = value;
            }
        }
    }

    public void saveTo(DataOutputStream out) throws IOException {
        out.writeLong(size);
        for (long i = 0; i < size; ++i) {
            out.writeLong(get(i));
        }
    }

    public static LLongArray loadFrom(DataInputStream in) throws IOException {
        long size = in.readLong();
        LLongArray array = new LLongArray(size);
        for (long i = 0; i < size; ++i) {
            array.set(i, in.readLong());
        }
        return array;
    }

}
