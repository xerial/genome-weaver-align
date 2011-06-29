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
// Client.java
// Since: 2011/04/25
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.parallel;

import org.msgpack.rpc.loop.EventLoop;
import org.utgenome.weaver.align.GenomeWeaverCommand;
import org.xerial.util.log.Logger;
import org.xerial.util.opt.Argument;
import org.xerial.util.opt.Option;

public class Client extends GenomeWeaverCommand
{
    private static Logger _logger = Logger.getLogger(Client.class);

    @Override
    public String name() {
        return "client";
    }

    @Override
    public String getOneLineDescription() {
        return "Launch a client command";
    }

    public static interface RPCInterface
    {
        String hello(String message);
    }

    @Option(symbol = "p", description = "port number. default=8990")
    private int    port   = 8990;

    @Option(symbol = "s", description = "server address. default=localhost")
    private String server = "localhost";

    @Argument()
    private String input  = "hello";

    @Override
    public void execute(String[] args) throws Exception {
        EventLoop loop = EventLoop.defaultEventLoop();

        org.msgpack.rpc.Client cli = new org.msgpack.rpc.Client("localhost", 8990, loop);
        RPCInterface iface = cli.proxy(RPCInterface.class);
        String message = iface.hello(input);

        _logger.info(String.format("Recieved a message: %s", message));
    }

}
