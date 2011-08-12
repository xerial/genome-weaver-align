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
// BidirectionalSuffixInterval.java
// Since: 2011/08/08
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

/**
 * Suffix array interval holder for bi-directional search
 * 
 * @author leo
 * 
 */
public class BidirectionalSuffixInterval implements SARange
{
    public final SuffixInterval forwardSi;
    public final SuffixInterval backwardSi;

    public BidirectionalSuffixInterval(SuffixInterval forwardSi, SuffixInterval backwardSi) {
        this.forwardSi = forwardSi;
        this.backwardSi = backwardSi;
    }

    @Override
    public SuffixInterval forwardSi() {
        return forwardSi;
    }

    @Override
    public SuffixInterval backwardSi() {
        return backwardSi;
    }

}
