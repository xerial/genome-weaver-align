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
// BinInGenome.java
// Since: 2011/05/31
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.db;

import java.util.Iterator;
import java.util.List;

import org.utgenome.gwt.utgb.client.bio.GenomeRange;
import org.utgenome.gwt.utgb.client.bio.Interval;

/**
 * Holds entries in an interval of a chromosome
 * 
 * @author leo
 * 
 */
public class BinInGenome<T extends GenomeRange> implements Iterable<T>
{
    public final String   chr;
    public final Interval range;
    private final List<T> data;

    public BinInGenome(String chr, Interval bin, List<T> data) {
        this.chr = chr;
        this.range = bin;
        this.data = data;
    }

    public void add(T elem) {
        this.data.add(elem);
    }

    public boolean overlaps(GenomeRange in) {
        return range.hasOverlap(in);
    }

    public List<T> data() {
        return data;
    }

    public Iterator<T> iterator() {
        return data.iterator();
    }

    @Override
    public String toString() {
        return String.format("%s:%s", chr, range);
    }
}
