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
// BWFiles.java
// Since: 2011/02/15
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.io.File;

/**
 * This class defines BWT related file names
 * 
 * @author leo
 * 
 */
public class BWTFiles
{
    private String prefix;

    public BWTFiles(String fastaFile) {
        if (fastaFile.endsWith("tar.gz"))
            this.prefix = fastaFile.substring(0, fastaFile.length() - "tar.gz".length() - 1);
        else
            this.prefix = fastaFile;
    }

    public File pacFileIndex() {
        return new File(prefix + ".i.silk");
    }

    public File iupacForward() {
        return new File(prefix + ".f.iupac");
    }

    public File iupacReverse() {
        return new File(prefix + ".r.iupac");
    }

    public File bwtForward() {
        return new File(prefix + ".f.bwt");
    }

    public File bwtReverse() {
        return new File(prefix + ".r.bwt");
    }

    public File sparseSuffixArrayForward() {
        return new File(prefix + ".f.ssa");
    }

    public File sparseSuffixArrayReverse() {
        return new File(prefix + ".r.ssa");
    }

    public File bwtWaveletForward() {
        return new File(prefix + ".f.bwt.wv");
    }

    public File bwtWaveletReverse() {
        return new File(prefix + ".r.bwt.wv");
    }

}
