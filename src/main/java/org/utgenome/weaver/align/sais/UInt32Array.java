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
package org.utgenome.weaver.align.sais;

import java.util.ArrayList;
import java.util.Iterator;

import org.utgenome.weaver.align.LSeq;

/**
 * Array capable to store at most 4G (4 x 1024 x 1024 x 1024) entries
 * 
 * @author leo
 * 
 */
public class UInt32Array implements LSeq, Iterable<Long>
{
    private static final int  B           = 30;            // bit length 
    private static final int  BLOCK_SIZE  = 1 << B;
    private static final long OFFSET_MASK = BLOCK_SIZE - 1;

    private final long        size;

    private ArrayList<int[]>  array;

    public UInt32Array(long size) {
        this.size = size;

        int container = (int) (size >>> B);
        int remainder = offset(size);

        array = new ArrayList<int[]>(container + 1);
        for (int i = 0; i < container; i++) {
            array.add(new int[BLOCK_SIZE]);
        }
        if (remainder > 0)
            array.add(new int[remainder]);
    }

    private int[] container(long index) {
        return array.get((int) (index >>> B));
    }

    private static int offset(long index) {
        return (int) (index & OFFSET_MASK);
    }

    public long textSize() {
        return size;
    }

    public long lookup(long index) {
        long v = container(index)[offset(index)] & 0xFFFFFFFFL;
        return v;
    }

    public void set(long index, long value) {
        container(index)[offset(index)] = (int) value;
    }

    @Override
    public long update(long index, long val) {
        long next = lookup(index) + val;
        set(index, next);
        return next;
    }

    @Override
    public String toString() {
        StringBuilder b = new StringBuilder();
        int i = 0;
        b.append("[");
        for (long each : this) {
            if (i++ > 0)
                b.append(", ");
            b.append(each);
        }
        b.append("]");
        return b.toString();
    }

    public long[] toArray() {
        long[] r = new long[(int) textSize()];
        for (long i = 0; i < textSize(); ++i) {
            r[(int) i] = lookup(i);
        }
        return r;
    }

    @Override
    public Iterator<Long> iterator() {
        return new Iterator<Long>() {

            private long cursor = 0;

            @Override
            public boolean hasNext() {
                return cursor < size;
            }

            @Override
            public Long next() {
                return lookup(cursor++);
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException("remove");
            }
        };
    }

}
