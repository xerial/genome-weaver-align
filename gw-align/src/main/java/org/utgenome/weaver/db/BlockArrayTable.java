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
// BlockArrayTable.java
// Since: 2011/04/06
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.db;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;

import org.utgenome.weaver.db.BlockArray.Block;
import org.xerial.db.DBException;
import org.xerial.db.sql.ResultSetHandler;
import org.xerial.db.sql.sqlite.SQLiteAccess;
import org.xerial.util.log.Logger;

/**
 * Manages {@link BlockArray} of chromosomes
 * 
 * @author leo
 * 
 */
public class BlockArrayTable
{
    private static Logger               _logger = Logger.getLogger(BlockArrayTable.class);

    private TreeMap<String, BlockArray> table   = new TreeMap<String, BlockArray>();

    public BlockArray getBlockArray(String chr) {
        if (!table.containsKey(chr))
            table.put(chr, new BlockArray());
        return table.get(chr);
    }

    public float get(String chr, int index) {
        return getBlockArray(chr).get(index);
    }

    public Set<String> keySet() {
        return table.keySet();
    }

    public void put(String chr, BlockArray block) {
        table.put(chr, block);
    }

    public void saveTo(File f) throws IOException {
        DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(f)));

        out.writeInt(table.size());
        for (String chr : table.keySet()) {
            BlockArray b = table.get(chr);
            byte[] key = chr.getBytes("UTF-8");

            out.writeInt(key.length);
            out.write(key);
            b.saveTo(out);
        }

        out.close();
    }

    public static BlockArrayTable loadFrom(File f) throws IOException {
        DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(f)));
        BlockArrayTable table = new BlockArrayTable();

        int numEntries = in.readInt();
        for (int i = 0; i < numEntries; ++i) {
            int chrNameByteSize = in.readInt();
            byte[] chrNameBytes = new byte[chrNameByteSize];
            in.read(chrNameBytes);
            String chr = new String(chrNameBytes, "UTF-8");
            BlockArray blockArray = BlockArray.loadFrom(in);
            table.put(chr, blockArray);
        }
        in.close();
        return table;
    }

    private static class ChromIndex
    {
        public int    trackId;
        public String chr;
    }

    public static BlockArrayTable loadFromSQLite(File f) throws IOException, DBException {
        SQLiteAccess db = new SQLiteAccess(f.getAbsolutePath());
        List<ChromIndex> chrSet = db.query("select track_id, value as chr from track where name=\"chrom\"",
                ChromIndex.class);

        BlockArrayTable t = new BlockArrayTable();
        for (ChromIndex eachChr : chrSet) {
            _logger.info("Loading " + eachChr.chr);
            final BlockArray blockArray = t.getBlockArray(eachChr.chr);
            db.query(String.format("select start, end, data_values from data where track_id = %d order by start",
                    eachChr.trackId), new ResultSetHandler<Void>() {
                @Override
                public Void handle(ResultSet rs) throws SQLException {
                    int start = rs.getInt("start");
                    int end = rs.getInt("end");
                    ObjectInputStream in = null;
                    try {
                        in = new ObjectInputStream(new GZIPInputStream(new ByteArrayInputStream(rs
                                .getBytes("data_values"))));
                        float[] data = (float[]) in.readObject();
                        blockArray.add(new Block(start, data));
                    }
                    catch (Exception e) {
                        e.printStackTrace();
                    }
                    finally {
                        if (in != null)
                            try {
                                in.close();
                            }
                            catch (IOException e) {
                                e.printStackTrace();
                            }
                    }
                    return null;
                }
            });
        }

        return t;

    }

}
