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
// PairedEndRead.java
// Since: 2011/04/28
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.record;

import org.xerial.util.StringUtil;

/**
 * Paired-end read sequence holder
 * 
 * @author leo
 * 
 */
public class PairedEndRead implements RawRead
{
    public final ReadSequence first;
    public final ReadSequence second;

    public PairedEndRead(ReadSequence first, ReadSequence second) {
        this.first = first;
        this.second = second;
    }

    @Override
    public String name() {
        return first.name;
    }

    @Override
    public String toString() {
        return StringUtil.join(new ReadSequence[] { first, second }, "\n");
    }

}
