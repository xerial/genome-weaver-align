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

import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

import org.xerial.util.StringUtil;
import org.xerial.util.log.Logger;
import org.xerial.util.opt.CommandHelpMessage;
import org.xerial.util.opt.CommandLauncher;

/**
 * Entry point of the Genome Weaver
 * 
 * @author leo
 * 
 */
public class GenomeWeaver
{
    private static Logger _logger = Logger.getLogger(GenomeWeaver.class);

    public static void execute(String arg) throws Exception {
        execute(StringUtil.tokenizeCommandLineArgument(arg));
    }

    public static void execute(String[] arg) throws Exception {
        CommandLauncher l = new CommandLauncher();
        CommandHelpMessage m = new CommandHelpMessage();
        m.defaultHeader = String.format("Genome Weaver: version %s", getVersion());
        l.setMessage(m);
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

    public static String getVersion() {
        String version = "(unknown)";
        try {
            // load the pom.xml file copied as a resource in utgb-core.jar
            String propertyName = "version";
            InputStream pomIn = GenomeWeaver.class
                    .getResourceAsStream("/META-INF/maven/org.utgenome.weaver/genome-weaver/pom.properties");
            try {
                if (pomIn == null) {
                    // If genome-pweaver is referenced in the workspace scope, use the
                    // genome-weaver/src/main/resources/genome-weaver.properties, which is created when genome-weaver is
                    // compiled
                    pomIn = GenomeWeaver.class.getResourceAsStream("/org/utgenome/weaver/genome-weaver.properties");
                    propertyName = "genome-weaver-version";
                }
                if (pomIn != null) {
                    Properties prop = new Properties();
                    prop.load(pomIn);
                    version = prop.getProperty(propertyName, version);
                }
            }
            finally {
                if (pomIn != null)
                    pomIn.close();
            }
        }
        catch (IOException e) {
            _logger.debug(e);
        }
        return version;
    }

}
