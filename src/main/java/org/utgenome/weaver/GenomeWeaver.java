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
// GenomeWeaver.java
// Since: 2011/02/10
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver;

import org.xerial.util.StringUtil;
import org.xerial.util.opt.CommandLauncher;

/**
 * Entry point of the Genome Weaver
 * 
 * @author leo
 * 
 */
public class GenomeWeaver
{
    public static void execute(String arg) throws Exception {
        execute(StringUtil.tokenizeCommandLineArgument(arg));
    }

    public static void execute(String[] arg) throws Exception {
        CommandLauncher l = new CommandLauncher();
        l.addCommandsIn(GenomeWeaver.class.getPackage(), true);
        l.execute(arg);
    }

    public static void main(String[] args) {
        try {
            GenomeWeaver.execute(args);
        }
        catch (Exception e) {
            e.printStackTrace(System.err);
        }
    }

}
