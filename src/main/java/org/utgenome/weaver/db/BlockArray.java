/*--------------------------------------------------------------------------
 *  Copyright 2008 utgenome.org
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
// tss-toolkit Project
//
// BlockArray.java
// Since: 2011/04/05
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.db;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.util.Map.Entry;
import java.util.TreeMap;

import org.utgenome.gwt.utgb.client.bio.Interval;
import org.xerial.snappy.Snappy;

/**
 * BlockArray stores a chain of blocks. Each block has a pair of (start pos,
 * data[]).
 * 
 * 
 * 
 * @author leo
 * 
 */
public class BlockArray
{

    public static class Block implements Comparable<Block>
    {
        public final int start;
        public float[]   data;

        public Block(int start) {
            this.start = start;
        }

        public Block(Interval bin) {
            this.start = bin.getStart();
            this.data = new float[bin.length()];
        }

        public Block(int start, float[] data) {
            this.start = start;
            this.data = data;
        }

        public int compareTo(Block o) {
            return this.start - o.start;
        }

        public void set(int index, float value) {
            int offset = index - start;
            if (offset < 0 || offset >= data.length) {
                throw new IndexOutOfBoundsException(String.format("index:%d", index));
            }
            this.data[offset] = value;
        }

        public float get(int index) {
            if ((index - start) < data.length)
                return data[index - start];
            else
                return 0;
        }

        public int getEnd() {
            return this.start + data.length;
        }
    }

    private TreeMap<Integer, Block> blockChain = new TreeMap<Integer, Block>();

    public BlockArray() {}

    /**
     * Return the element at the specified index.
     * 
     * @param index
     *            index in the array
     * @return element value. 0 will be returned if no entry is found at the
     *         index.
     */
    public float get(int index) {
        Entry<Integer, Block> floorEntry = blockChain.floorEntry(index);
        if (floorEntry == null)
            return 0;

        return floorEntry.getValue().get(index);
    }

    public void add(Block block) {
        blockChain.put(block.start, block);
    }

    public int getMinPosition() {
        Entry<Integer, Block> firstEntry = blockChain.firstEntry();
        if (firstEntry == null)
            return 0;
        else {

            return firstEntry.getValue().start;
        }
    }

    public int getMaxLength() {
        Entry<Integer, Block> lastEntry = blockChain.lastEntry();
        if (lastEntry == null)
            return 0;

        return lastEntry.getValue().getEnd();
    }

    public void saveTo(DataOutputStream out) throws IOException {

        out.writeInt(blockChain.size());
        for (int start : blockChain.keySet()) {
            Block block = blockChain.get(start);
            ByteArrayOutputStream buf = new ByteArrayOutputStream();
            DataOutputStream dout = new DataOutputStream(buf);
            for (float f : block.data)
                dout.writeFloat(f);
            dout.close();
            byte[] compressed = Snappy.compress(buf.toByteArray());

            // Write an entry (start, compressed length, compresesd data...)
            out.writeInt(block.start);
            out.writeInt(compressed.length);
            out.write(compressed);
        }

    }

    public static BlockArray loadFrom(DataInputStream in) throws IOException {

        BlockArray array = new BlockArray();
        try {
            int numBlocks = in.readInt();
            for (int k = 0; k < numBlocks; k++) {
                int start = in.readInt();
                int compressedDataLength = in.readInt();
                byte[] compressedFloatData = new byte[compressedDataLength];
                in.read(compressedFloatData);

                byte[] uncompressedFloat = Snappy.uncompress(compressedFloatData);
                DataInputStream floatIn = new DataInputStream(new ByteArrayInputStream(uncompressedFloat));
                float[] data = new float[uncompressedFloat.length / 4];
                for (int i = 0; i < data.length; ++i)
                    data[i] = floatIn.readFloat();

                Block block = new Block(start);
                block.data = data;
                array.add(block);
            }
            return array;
        }
        catch (EOFException e) {
            return array;
        }

    }
}
