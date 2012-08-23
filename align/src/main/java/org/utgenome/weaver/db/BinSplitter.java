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
// BinSplitter.java
// Since: 2011/05/27
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.db;

import java.util.ArrayList;

import org.utgenome.gwt.utgb.client.bio.ChrInterval;
import org.utgenome.gwt.utgb.client.bio.GenomeRange;
import org.utgenome.gwt.utgb.client.bio.GenomeRangeDataSet;
import org.utgenome.gwt.utgb.client.bio.Interval;
import org.xerial.util.ObjectHandler;

/**
 * Receives a sorted input of {@link GenomeRange} data, then create bins of
 * {@link GenomeRangeDataSet}
 * 
 * @author leo
 * 
 */
public class BinSplitter<T extends ChrInterval> implements ObjectHandler<T>
{

    private final int                           binSize;
    private final ObjectHandler<BinInGenome<T>> out;

    private int                                 cursor;

    private BinInGenome<T>                      currentBin;

    public BinSplitter(int binSize, ObjectHandler<BinInGenome<T>> out) {
        this.binSize = binSize;
        this.out = out;
    }

    public void init() throws Exception {
        out.init();
    }

    public void handle(T input) throws Exception {
        int start = input.getStart();
        int binIndex = start / binSize;
        if (currentBin == null) {
            int binStart = (start / binSize) * binSize;
            currentBin = new BinInGenome<T>(input.chr, new Interval(binStart, binStart + binSize), new ArrayList<T>());
        }

        boolean isInSameChr = currentBin.chr.equals(input.chr);
        if (isInSameChr && currentBin.overlaps(input)) {
            currentBin.add(input);
        }
        else {
            // sweep the current bin
            out.handle(currentBin);

            // Create a new bin
            int binStart = (start / binSize) * binSize;
            BinInGenome<T> newBin = new BinInGenome<T>(input.chr, new Interval(binStart, binStart + binSize),
                    new ArrayList<T>());

            if (isInSameChr) {
                // Put the entries which also overlap with the new bin
                for (T each : currentBin) {
                    if (newBin.overlaps(each)) {
                        newBin.add(each);
                    }
                }
            }
            newBin.add(input);
            currentBin = newBin;
        }
    }

    public void finish() throws Exception {
        if (currentBin != null)
            out.handle(currentBin);

        out.finish();
    }

}
