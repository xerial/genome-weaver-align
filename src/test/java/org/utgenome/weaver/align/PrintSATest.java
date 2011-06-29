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
// PrintSATest.java
// Since: 2011/02/17
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import org.junit.Test;
import org.utgenome.weaver.GenomeWeaver;

public class PrintSATest
{
    @Test
    public void print() throws Exception {
        GenomeWeaver.execute("sa-debug ACGTCCTTATACTTTATCTCATGGGATA");
    }

    @Test
    public void reverse() throws Exception {
        GenomeWeaver.execute("sa-debug -r ACGTCCTTATACTTTATCTCATGGGATA");
    }

    @Test
    public void sample2() throws Exception {
        GenomeWeaver.execute("sa-debug AACCTATCTATACCCCGGGGAATTATATACGTCCTTATACTTTATCTCATGGGATA");
    }

    @Test
    public void sample3() throws Exception {
        GenomeWeaver.execute("sa-debug AGGAGC");
    }

    @Test
    public void sample4() throws Exception {
        GenomeWeaver.execute("sa-debug AGGAGC -r");
    }

}
