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
// RawRead.java
// Since: 2011/04/28
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.record;

import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.SOLiDColorSequence;

public interface Read
{
    public String name();

    public int getNumReadFragment();

    public ACGTSequence getRead(int index);

    public SOLiDColorSequence getColorRead(int index);

    public boolean isLetterSpace();

    public boolean isColorSpace();

    public String getQual(int index);
}
