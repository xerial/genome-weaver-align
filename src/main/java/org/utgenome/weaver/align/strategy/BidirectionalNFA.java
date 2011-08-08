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
// BidirectionalNFA.java
// Since: Aug 5, 2011
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.strategy;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.PriorityQueue;

import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.SuffixInterval;
import org.utgenome.weaver.parallel.Reporter;
import org.xerial.lens.SilkLens;
import org.xerial.util.Optional;
import org.xerial.util.log.Logger;

/**
 * NFA for forward and backward traversal of suffix arrays
 * 
 * @author leo
 * 
 */
public class BidirectionalNFA
{
    private static Logger         _logger              = Logger.getLogger(BidirectionalNFA.class);

    private final FMIndexOnGenome fmIndex;
    private final ACGTSequence    qF;
    private final ACGTSequence    qC;
    private final int             m;                                                              // query length
    private final int             numAllowedMismatches = 1;

    public BidirectionalNFA(FMIndexOnGenome fmIndex, ACGTSequence query) {
        this.fmIndex = fmIndex;
        this.qF = query;
        this.qC = query.complement();
        this.m = (int) query.textSize();
    }

    private static class CursorContainer
    {
        private List<Optional<List<BidirectionalCursor>>> containerF;
        private List<Optional<List<BidirectionalCursor>>> containerR;

        public CursorContainer(int m) {
            containerF = new ArrayList<Optional<List<BidirectionalCursor>>>(m);
            containerR = new ArrayList<Optional<List<BidirectionalCursor>>>(m);
            for (int i = 0; i < m; ++i) {
                containerF.add(new Optional<List<BidirectionalCursor>>());
                containerR.add(new Optional<List<BidirectionalCursor>>());
            }
        }

        private List<Optional<List<BidirectionalCursor>>> getContainer(Strand strand) {
            return strand.isForward() ? containerF : containerR;
        }

        public void add(BidirectionalCursor c) {
            Optional<List<BidirectionalCursor>> l = getContainer(c.strand).get(c.getIndex());
            if (l.isUndefined()) {
                l.set(new ArrayList<BidirectionalCursor>());
            }
            l.get().add(c);
        }

        @Override
        public String toString() {
            return String.format("%s\n%s", containerF, containerR);
        }
    }

    public static class ExactMatch
    {
        public final Strand         strand;
        public final SuffixInterval si;
        public final int            numMismatches;

        public ExactMatch(Strand strand, SuffixInterval si, int numMismatches) {
            this.strand = strand;
            this.si = si;
            this.numMismatches = numMismatches;
        }
    }

    public void align() throws Exception {
        bidirectionalAlign(new Reporter() {
            @Override
            public void emit(Object result) throws Exception {
                _logger.debug(SilkLens.toSilk(result));
            }
        });
    }

    /**
     * @param out
     * @throws Exception
     */
    public void bidirectionalAlign(Reporter out) throws Exception {

        // Prefer the cursor with smaller remaining bases
        PriorityQueue<BidirectionalCursor> queue = new PriorityQueue<BidirectionalCursor>(11,
                new Comparator<BidirectionalCursor>() {
                    @Override
                    public int compare(BidirectionalCursor o1, BidirectionalCursor o2) {
                        return o1.getRemainingBases() - o2.getRemainingBases();
                    }
                });

        // Row-wise simulation of NFA
        CursorContainer state = new CursorContainer(m + 1);

        BidirectionalCursor initF = new BidirectionalCursor(qF, Strand.FORWARD, SearchDirection.Forward,
                fmIndex.wholeSARange(), null, 0, 0);
        BidirectionalCursor initC = new BidirectionalCursor(qC, Strand.REVERSE, SearchDirection.Forward,
                fmIndex.wholeSARange(), null, 0, 0);
        queue.add(initF);
        queue.add(initC);
        state.add(initF);
        state.add(initC);

        for (int k = 0; k < numAllowedMismatches; ++k) {

            while (!queue.isEmpty()) {
                BidirectionalCursor c = queue.poll();
                if (c.getRemainingBases() == 0) {
                    // Found a match
                    out.emit(c);
                    continue;
                }

                BidirectionalCursor next = c.next(fmIndex);
                if (next == null)
                    continue; // no match
                queue.add(next);

                state.add(next);
            }
        }

    }

}
