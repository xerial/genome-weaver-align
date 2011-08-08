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

import org.utgenome.weaver.align.ACGT;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.AlignmentScoreConfig;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.Strand;
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
    private static Logger              _logger   = Logger.getLogger(BidirectionalNFA.class);

    private final FMIndexOnGenome      fmIndex;
    private final ACGTSequence         qF;
    private final ACGTSequence         qC;
    private final int                  m;                                                   // query length
    private final AlignmentScoreConfig config;
    private final Reporter             out;
    private int                        bestScore = -1;

    public BidirectionalNFA(FMIndexOnGenome fmIndex, ACGTSequence query) {
        this(fmIndex, query, new Reporter() {
            @Override
            public void emit(Object result) throws Exception {
                _logger.debug(SilkLens.toSilk(result));
            }
        });
    }

    public BidirectionalNFA(FMIndexOnGenome fmIndex, ACGTSequence query, Reporter out) {
        this.fmIndex = fmIndex;
        this.qF = query;
        this.qC = query.complement();
        this.m = (int) query.textSize();
        this.config = new AlignmentScoreConfig();
        this.out = out;
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

    private static class CursorQueue extends PriorityQueue<BidirectionalCursor>
    {
        /**
         * 
         */
        private static final long serialVersionUID = 1L;

        public CursorQueue() {
            super(11, new Comparator<BidirectionalCursor>() {
                @Override
                public int compare(BidirectionalCursor o1, BidirectionalCursor o2) {
                    return o1.getRemainingBases() - o2.getRemainingBases();
                }
            });
        }

    }

    /**
     * @param out
     * @throws Exception
     */
    public void align() throws Exception {

        // Prefer the cursor with smaller remaining bases
        CursorQueue queue = new CursorQueue();
        CursorQueue nextQueue = new CursorQueue();

        // Row-wise simulation of NFA
        CursorContainer state = new CursorContainer(m + 1);

        BidirectionalCursor initF = new BidirectionalCursor(Score.initial(), qF, Strand.FORWARD,
                SearchDirection.Forward, fmIndex.wholeSARange(), null, 0, 0);
        BidirectionalCursor initC = new BidirectionalCursor(Score.initial(), qC, Strand.REVERSE,
                SearchDirection.Forward, fmIndex.wholeSARange(), null, 0, 0);
        queue.add(initF);
        queue.add(initC);
        state.add(initF);
        state.add(initC);

        // k=0;
        while (!queue.isEmpty()) {
            BidirectionalCursor c = queue.poll();
            if (!continueSearch(c))
                continue;

            // Proceed to next base
            BidirectionalCursor next = c.next(fmIndex, config);
            if (next != null) {
                queue.add(next);
                state.add(next);
            }
            else {
                // pass to next layer
                nextQueue.add(c);
            }

        }

        // k >= 1
        for (int k = 1; k < config.maximumEditDistances; ++k) {
            // Transit the states to the next row
            _logger.debug("transit k from %d to %d", k - 1, k);
            queue = nextQueue;
            nextQueue = new CursorQueue();

            while (!queue.isEmpty()) {
                BidirectionalCursor c = queue.poll();
                if (!continueSearch(c))
                    continue;

                if (c.score.layer() == k) {
                    // search for exact match
                    BidirectionalCursor next = c.next(fmIndex, config);
                    if (next != null)
                        queue.add(next);
                }
                else {
                    // extend the search with mismatches
                    ACGT nextBase = c.nextACGT();
                    for (ACGT ch : ACGT.exceptN) {
                        if (nextBase == ch)
                            continue;
                        BidirectionalCursor next = c.next(fmIndex, ch, config);
                        if (next != null)
                            queue.add(next);
                    }
                }
            }

        }

    }

    public boolean continueSearch(BidirectionalCursor c) throws Exception {
        if (c.getUpperBoundOfScore(config) < bestScore) {
            // This cursor never produces a better alignment. Skip 
            return false;
        }

        if (c.getRemainingBases() == 0) {
            // Update the best score
            if (c.score.score > bestScore) {
                bestScore = c.score.score;
            }
            // Found a match
            out.emit(c);
            return false;
        }
        return true;
    }

}
