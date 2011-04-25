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
// Server.java
// Since: 2011/04/25
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.parallel;

import org.msgpack.rpc.loop.EventLoop;
import org.utgenome.weaver.align.GenomeWeaverCommand;
import org.xerial.util.log.Logger;
import org.xerial.util.opt.Option;

public class Server extends GenomeWeaverCommand
{
    private static Logger _logger = Logger.getLogger(Server.class);

    @Override
    public String name() {
        return "server";
    }

    @Override
    public String getOneLineDescription() {
        return "Launch a genome-weaver server";
    }

    public class JobManager extends org.msgpack.rpc.Server
    {
        public JobManager(EventLoop loop) {
            super(loop);
        }

        public String hello(String message) {
            return String.format("Hello %s!", message);
        }

    }

    @Option(symbol = "p", description = "listen port. default = 8990")
    private int port = 8990;

    @Override
    public void execute(String[] args) throws Exception {

        EventLoop loop = EventLoop.defaultEventLoop();
        JobManager jobManager = new JobManager(loop);
        jobManager.serve(jobManager);
        jobManager.listen(port);

        loop.join();
    }

}
