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
// LaunchJava.java
// Since: 2011/02/16
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.net.URL;
import java.util.ArrayList;
import java.util.Map.Entry;
import java.util.Properties;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.xerial.util.StringUtil;
import org.xerial.util.log.Logger;
import org.xerial.util.opt.CommandModule;
import org.xerial.util.opt.GlobalCommandOption;

public class LaunchJava implements CommandModule
{
    private static Logger _logger = Logger.getLogger(LaunchJava.class);

    @Override
    public String name() {
        return "e";
    }

    @Override
    public String getOneLineDescription() {
        return "launch a command in a child process";
    }

    @Override
    public Object getOptionHolder() {
        return this;
    }

    @Override
    public void execute(String[] args) throws Exception {
        if (args.length <= 0) {
            _logger.error("Empty command");
            return;
        }

        CommandExecutor.exec(concatenate(args), CommandExecutor.prepareEnvironmentVariables(), null);
    }

    private static String concatenate(String[] args) {
        ArrayList<String> cmd = new ArrayList<String>();
        for (String each : args) {
            if (each.contains(" ")) {
                each = StringUtil.doubleQuote(each);
            }
            cmd.add(each);
        }
        return StringUtil.join(cmd, " ");
    }

    private static abstract class ProcessOutputReader implements Runnable
    {
        private final BufferedReader reader;

        public ProcessOutputReader(InputStream in) {
            reader = new BufferedReader(new InputStreamReader(in));
        }

        public abstract void output(String line);

        public void run() {
            try {
                String line;
                while ((line = reader.readLine()) != null) {
                    output(line);
                }
            }
            catch (IOException e) {
                // If the process is already terminated, IOException (bad file descriptor) might be reported.
                _logger.debug(e);
            }
        }
    }

    private static class ProcessInputWriter implements Runnable
    {
        private final PrintWriter    writer;
        private final BufferedReader systemIn;

        public ProcessInputWriter(OutputStream out) {
            writer = new PrintWriter(out);
            systemIn = new BufferedReader(new InputStreamReader(System.in));
        }

        @Override
        public void run() {
            try {
                String line;
                while ((line = systemIn.readLine()) != null) {
                    writer.println(line);
                }
            }
            catch (IOException e) {
                _logger.debug(e);
            }

        }

    }

    public static class CommandExecutor
    {
        final ExecutorService threadManager = Executors.newFixedThreadPool(3);
        Process               proc          = null;
        Future< ? >           stdoutReader;
        Future< ? >           stderrReader;

        //Future< ? >           stdInWriter;

        private void dispose() {
            if (proc != null) {
                proc.destroy();
                proc = null;
            }

            threadManager.shutdown();
            try {
                while (!threadManager.awaitTermination(1L, TimeUnit.SECONDS)) {}
            }
            catch (InterruptedException e) {
                _logger.error(e);
            }
        }

        public int execCommand(String commandLine, String[] envp, File workingDir) throws IOException {
            try {
                if (_logger.isDebugEnabled())
                    _logger.debug(commandLine);

                proc = Runtime.getRuntime().exec(commandLine, envp, workingDir);

                // pipe the program's stdout and stderr to the logger
                stdoutReader = threadManager.submit(new ProcessOutputReader(proc.getInputStream()) {
                    @Override
                    public void output(String line) {
                        System.out.println(line);
                    }
                });
                stderrReader = threadManager.submit(new ProcessOutputReader(proc.getErrorStream()) {
                    @Override
                    public void output(String line) {
                        System.err.println(line);
                    }
                });

                //stdInWriter = threadManager.submit(new ProcessInputWriter(proc.getOutputStream()));

                int ret = proc.waitFor();
                return ret;
            }
            catch (InterruptedException e) {
                _logger.error(e);
                return 0;
            }
            finally {
                dispose();
            }

        }

        public static int exec(String commandLine) throws IOException {
            return exec(commandLine, null, null);
        }

        public static int exec(String commandLine, String[] envp, File workingDir) throws IOException {
            CommandExecutor e = new CommandExecutor();
            return e.execCommand(commandLine, envp, workingDir);
        }

        public static String[] prepareEnvironmentVariables() {
            Properties env = new Properties();
            for (Entry<String, String> eachEnv : System.getenv().entrySet()) {
                env.setProperty(eachEnv.getKey(), eachEnv.getValue());
            }
            if (!env.contains("JAVA_HOME") || env.getProperty("JAVA_HOME").contains("jre")) {
                env.setProperty("JAVA_HOME", System.getProperty("java.home"));
            }

            String[] envp = new String[env.size()];
            int index = 0;
            for (Object each : env.keySet()) {
                String key = each.toString();
                envp[index++] = String.format("%s=%s", key, env.getProperty(key));
            }

            _logger.trace("environment variables: " + env);
            return envp;
        }

    }

    @Override
    public void printUsage() throws Exception {
        System.out.println("[usage]");
        System.out.println("$ genome-weaver e (command line...)");
    }

    @Override
    public URL getHelpMessageResource() {
        return null;
    }

    @Override
    public void execute(GlobalCommandOption globalOption, String[] args) throws Exception {
        execute(args);
    }

}
