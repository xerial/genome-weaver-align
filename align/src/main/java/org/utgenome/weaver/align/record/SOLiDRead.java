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
// SOLiDRead.java
// Since: 2011/10/18
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.record;

import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.SOLiDColorSequence;

public class SOLiDRead implements Read
{
    public final String             name;
    public final SOLiDColorSequence read;
    public final String             qual;

    public SOLiDRead(String name, SOLiDColorSequence read, String qual) {
        this.name = name;
        this.read = read;
        this.qual = qual;
    }

    @Override
    public String name() {
        return name;
    }

    @Override
    public int getNumReadFragment() {
        return 1;
    }

    @Override
    public ACGTSequence getRead(int index) {
        return null;
    }

    @Override
    public SOLiDColorSequence getColorRead(int index) {
        if (index > 1)
            throw new IndexOutOfBoundsException("invalid index " + index);

        return read;
    }

    @Override
    public boolean isLetterSpace() {
        return false;
    }

    @Override
    public boolean isColorSpace() {
        return true;
    }

    @Override
    public String getQual(int index) {
        return qual;
    }

}
