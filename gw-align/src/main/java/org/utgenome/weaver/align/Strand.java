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
// Strand.java
// Since: 2011/02/16
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

public enum Strand {
    FORWARD("+", 0), REVERSE("-", 1);

    public final String symbol;
    public final int    index;

    private Strand(String symbol, int index) {
        this.symbol = symbol;
        this.index = index;
    }

    public static Strand decode(int index) {
        return index == 0 ? FORWARD : REVERSE;
    }

    public static Strand toStrand(char plusOrMinus) {
        if (plusOrMinus == '+')
            return FORWARD;
        else
            return REVERSE;
    }

    public boolean isForward() {
        return this == FORWARD;
    }

    public Strand opposite() {
        if (this == FORWARD)
            return REVERSE;
        else
            return FORWARD;
    }

}
