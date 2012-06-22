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
// GeneSet.java
// Since: 2011/01/13
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.db;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;

import org.utgenome.UTGBException;
import org.utgenome.format.bed.BED2SilkReader;
import org.utgenome.gwt.utgb.client.bio.BEDGene;
import org.utgenome.gwt.utgb.client.bio.GenomeRange;
import org.utgenome.gwt.utgb.client.bio.Interval;
import org.utgenome.gwt.utgb.client.canvas.IntervalTree;
import org.xerial.lens.SilkLens;
import org.xerial.util.ArrayDeque;
import org.xerial.util.ObjectHandler;
import org.xerial.util.log.Logger;

/**
 * Priority search tree for genome sequences, consisting of several chromosomes.
 * 
 * @author leo
 * 
 */
public class PrioritySearchTreeInGenome<T extends GenomeRange> implements Iterable<T>
{

    private static Logger                    _logger = Logger.getLogger(PrioritySearchTreeInGenome.class);

    private HashMap<String, IntervalTree<T>> table   = new HashMap<String, IntervalTree<T>>();

    public Set<String> getChrSet() {
        return table.keySet();
    }

    private int numEntries = 0;

    public int numEntries() {
        return numEntries;
    }

    public IntervalTree<T> get(String chr) {
        if (!table.containsKey(chr))
            table.put(chr, new IntervalTree<T>());

        return table.get(chr);
    }

    public void put(String chr, T gene) {
        IntervalTree<T> intervalTree = get(chr);
        intervalTree.add(gene);
        numEntries++;
    }

    public List<T> overlapQuery(String chr, int pos, Strand strand) {
        IntervalTree<T> t = get(chr);
        List<T> result = new ArrayList<T>();
        for (T each : t.overlapQuery(pos)) {
            if (strand == Strand.PLUS) {
                if (each.isSense())
                    result.add(each);
            }
            else {
                if (!each.isSense())
                    result.add(each);
            }

        }
        return result;
    }

    public List<T> overlapQuery(String chr, Interval queryRange, Strand strand) {
        IntervalTree<T> t = get(chr);
        List<T> result = new ArrayList<T>();
        for (T each : t.overlapQuery(queryRange)) {
            if (strand == Strand.PLUS) {
                if (each.isSense())
                    result.add(each);
            }
            else {
                if (!each.isSense())
                    result.add(each);
            }

        }
        return result;
    }

    public T nearestNeighbour(String chr, T base, Strand strand, int searchRange) {
        int pos = base.getStart();
        Interval queryRange = new Interval(pos - searchRange, pos + searchRange);
        T nearestEntry = null;
        for (T each : overlapQuery(chr, queryRange, strand)) {
            if (each == base)
                continue;

            if (each.getStart() == pos)
                continue;

            if (nearestEntry == null) {
                nearestEntry = each;
                continue;
            }

            int maxDist = Math.abs(nearestEntry.getStart() - pos);
            int currentDist = Math.abs(each.getStart() - pos);
            if (currentDist < maxDist) {
                nearestEntry = each;
            }
        }

        return nearestEntry;
    }

    public void loadBED(final String file) throws UTGBException, IOException {
        Reader r = null;
        try {
            r = new BED2SilkReader(new BufferedReader(new FileReader(file)));
            SilkLens.findFromSilk(new BED2SilkReader(new BufferedReader(new FileReader(file))), "gene", BEDGene.class,
                    new ObjectHandler<BEDGene>() {
                        int numEntries = 0;

                        public void init() throws Exception {
                            _logger.info("loading " + file);
                        }

                        @SuppressWarnings("unchecked")
                        public void handle(BEDGene input) throws Exception {
                            put(input.chr, (T) input);
                            numEntries++;
                        }

                        public void finish() throws Exception {
                            _logger.info(String.format("loaded %d entries", numEntries));
                        }
                    });
        }
        catch (Exception e) {
            throw UTGBException.convert(e);
        }
        finally {
            if (r != null)
                r.close();
        }

    }

    public Iterator<T> iterator() {
        return new Iterator<T>() {

            Iterator<String> chrCursor  = getChrSet().iterator();
            Iterator<T>      genesInChr = null;

            ArrayDeque<T>    queue      = new ArrayDeque<T>();

            public boolean hasNext() {

                if (!queue.isEmpty()) {
                    return true;
                }

                if (genesInChr == null) {
                    if (!chrCursor.hasNext())
                        return false;

                    String chr = chrCursor.next();

                    IntervalTree<T> intervalTree = get(chr);
                    genesInChr = intervalTree.iterator();
                }

                if (genesInChr.hasNext()) {
                    queue.add(genesInChr.next());
                    return true;
                }
                else {
                    genesInChr = null;
                    return hasNext();
                }
            }

            public T next() {
                if (!hasNext()) {
                    throw new NoSuchElementException();
                }

                return queue.pollFirst();
            }

            public void remove() {

            }

        };
    }
}
