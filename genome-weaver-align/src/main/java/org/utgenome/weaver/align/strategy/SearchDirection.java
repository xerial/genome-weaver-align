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
// SearchDirection.java
// Since: 2011/08/08
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.strategy;

/**
 * Search direction of a read. For a read [0, m), forward and
 * bidirectional-forward search traces the read from 0 to m, and backward search
 * traces the read from m-1 to 0.
 * 
 * @author leo
 * 
 */
public enum SearchDirection {
    Forward("F", 0, true), Backward("B", 1, false), BidirectionalForward("BF", 2, true);

    private final static SearchDirection[] lookupTable = { Forward, Backward, BidirectionalForward, Forward };

    public final String                    symbol;
    public final boolean                   isForward;
    public final int                       index;

    private SearchDirection(String symbol, int index, boolean isForward) {
        this.symbol = symbol;
        this.index = index;
        this.isForward = isForward;
    }

    public static SearchDirection decode(int index) {
        return lookupTable[index & 0x03];
    }

}
