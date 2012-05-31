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
// MismatchPosition.java
// Since: 2011/02/17
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.util.Arrays;

public class MismatchPosition
{
    private final static int B = 32;
    private final int[]      bitVector;
    private final int        size;

    public MismatchPosition(int size) {
        this.size = size;
        int blockSize = (size + B - 1) / B;
        bitVector = new int[blockSize];
        Arrays.fill(bitVector, 0);
    }

    private MismatchPosition(int[] bitVector, int size) {
        this.bitVector = bitVector;
        this.size = size;
    }

    public boolean get(int index) {
        int block = index / B;
        int offset = index % B;
        return (bitVector[block] & (0x80000000 >>> offset)) != 0;
    }

    public static MismatchPosition oneBitInstance(int size, int index) {
        int blockSize = (size + B - 1) / B;
        int[] bitVector = new int[blockSize];
        Arrays.fill(bitVector, 0);

        int block = index / B;
        int offset = index % B;
        bitVector[block] |= 0x80000000 >>> offset;
        return new MismatchPosition(bitVector, size);
    }

    public MismatchPosition copyAndSet(int index) {
        int[] newBitVector = new int[bitVector.length];
        System.arraycopy(bitVector, 0, newBitVector, 0, bitVector.length);
        int block = index / B;
        int offset = index % B;
        newBitVector[block] |= 0x80000000 >>> offset;
        return new MismatchPosition(newBitVector, size);
    }

    public MismatchPosition copyAndReset(int index) {
        int[] newBitVector = new int[bitVector.length];
        System.arraycopy(bitVector, 0, newBitVector, 0, bitVector.length);
        int block = index / B;
        int offset = index % B;
        newBitVector[block] &= ~(0x80000000 >>> offset);
        return new MismatchPosition(newBitVector, size);
    }

    @Override
    public String toString() {
        StringBuilder s = new StringBuilder();
        for (int i = 0; i < size; ++i) {
            s.append(get(i) ? "1" : "0");
        }
        return s.toString();
    }
}
