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

import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.SOLiDColorSequence;
import org.xerial.util.StringUtil;

/**
 * Paired-end read sequence holder
 * 
 * @author leo
 * 
 */
public class PairedEndRead implements Read
{
    public final Read first;
    public final Read second;

    public PairedEndRead(Read first, Read second) {
        this.first = first;
        this.second = second;
    }

    @Override
    public String name() {
        return first.name();
    }

    @Override
    public String toString() {
        return StringUtil.join(new Read[] { first, second }, "\n");
    }

    @Override
    public int getNumReadFragment() {
        return 1;
    }

    @Override
    public ACGTSequence getRead(int index) {
        switch (index) {
        case 0:
            return first.getRead(0);
        case 1:
            return second.getRead(0);
        default:
            throw new ArrayIndexOutOfBoundsException(index);
        }
    }

    @Override
    public String getQual(int index) {
        switch (index) {
        case 0:
            return first.getQual(0);
        case 1:
            return second.getQual(0);
        default:
            throw new ArrayIndexOutOfBoundsException(index);
        }
    }

    @Override
    public boolean isLetterSpace() {
        return true;
    }

    @Override
    public boolean isColorSpace() {
        return false;
    }

    @Override
    public SOLiDColorSequence getColorRead(int index) {
        return null;
    }
}
