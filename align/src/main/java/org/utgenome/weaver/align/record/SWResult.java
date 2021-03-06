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
// SWResult.java
// Since: 2011/08/22
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.record;

import org.xerial.lens.JSONLens;

/**
 * bit-parallel smith-waterman alignment result
 * 
 * @author leo
 * 
 */
public class SWResult
{
    public final int tailPos;
    public final int diff;

    public SWResult(int tailPos, int diff) {
        this.tailPos = tailPos;
        this.diff = diff;
    }

    @Override
    public String toString() {
        return JSONLens.toJSON(this);
    }
}
