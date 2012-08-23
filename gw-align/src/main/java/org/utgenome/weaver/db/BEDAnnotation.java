/*--------------------------------------------------------------------------
 *  Copyright 2009 utgenome.org
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
// utgb-core Project
//
// BEDGene.java
// Since: 2010/09/02
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.db;

import org.utgenome.gwt.utgb.client.bio.ChrInterval;

/**
 * Representing each gene line of BED format
 * 
 * @author yoshimura
 * 
 */
public class BEDAnnotation extends ChrInterval
{
    private static final long serialVersionUID = 1L;

    public String             name;
    public float              score            = 0;
    public String             strand;

    @Override
    public boolean isSense() {
        return "+".equals(strand);
    }

    @Override
    public boolean isAntiSense() {
        return !isSense();
    }

}
