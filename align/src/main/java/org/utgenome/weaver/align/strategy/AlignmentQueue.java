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
// AlignmentQueue.java
// Since: 2011/08/08
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.strategy;

import java.util.Comparator;
import java.util.PriorityQueue;

import org.utgenome.weaver.align.AlignmentScoreConfig;

public class AlignmentQueue
{
    final PriorityQueue<BWAState>     queue;
    private final AlignmentScoreConfig config;
    int                                bestScore;
    int                                pushCount = 0;

    public AlignmentQueue(AlignmentScoreConfig config) {
        this.config = config;
        this.queue = new PriorityQueue<BWAState>(11, new Comparator<BWAState>() {
            @Override
            public int compare(BWAState o1, BWAState o2) {
                // If the upper bound of the score is larger than the other, search it first
                int diff = o2.getUpperBoundOfScore(AlignmentQueue.this.config)
                        - o1.getUpperBoundOfScore(AlignmentQueue.this.config);
                if (diff != 0)
                    return diff;

                return o1.getRemainingBases() - o2.getRemainingBases();
            }
        });
    }

    @Override
    public String toString() {
        return String.format("size:%d, best:%d", queue.size(), bestScore);
    }

    public BWAState poll() {
        return queue.poll();
    }

    public boolean isEmpty() {
        return queue.isEmpty();
    }

    public boolean add(BWAState e) {
        if (e.getUpperBoundOfScore(config) < bestScore)
            return false;

        if (e.isFinished() && e.score.score > bestScore)
            bestScore = e.score.score;

        pushCount++;
        return queue.add(e);
    }

}
