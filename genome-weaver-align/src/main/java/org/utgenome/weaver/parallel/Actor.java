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
// Actor.java
// Since: 2011/10/19
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.parallel;

import java.lang.reflect.Constructor;
import java.lang.reflect.Method;
import java.net.InetSocketAddress;
import java.util.HashMap;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.jboss.netty.bootstrap.ServerBootstrap;
import org.jboss.netty.channel.Channel;
import org.jboss.netty.channel.ChannelFuture;
import org.jboss.netty.channel.ChannelHandlerContext;
import org.jboss.netty.channel.ChannelPipeline;
import org.jboss.netty.channel.ChannelPipelineFactory;
import org.jboss.netty.channel.Channels;
import org.jboss.netty.channel.MessageEvent;
import org.jboss.netty.channel.SimpleChannelUpstreamHandler;
import org.jboss.netty.channel.socket.nio.NioServerSocketChannelFactory;
import org.xerial.util.log.Logger;

public class Actor
{
    public static class RemoteServer
    {

        private static Logger               _logger    = Logger.getLogger(Actor.RemoteServer.class);

        public final String                 hostname;
        public final int                    port;
        private boolean                     isStarted  = false;
        private HashMap<String, Class< ? >> actorTable = new HashMap<String, Class< ? >>();

        public RemoteServer(String hostname, int port) {
            this.hostname = hostname;
            this.port = port;
        }

        public void start() {
            ExecutorService serverThread = Executors.newSingleThreadExecutor();
            Future<Boolean> ret = serverThread.submit(new Callable<Boolean>() {
                @Override
                public Boolean call() throws Exception {
                    ChannelFuture server = launchServer();
                    isStarted = true;
                    server.awaitUninterruptibly();
                    return true;
                }
            });
        }

        private ChannelFuture launchServer() {
            _logger.info("Start up a server %s:%s", hostname, port);

            ServerBootstrap bootstrap = new ServerBootstrap(new NioServerSocketChannelFactory(
                    Executors.newCachedThreadPool(), Executors.newCachedThreadPool()));

            bootstrap.setPipelineFactory(new ChannelPipelineFactory() {
                @Override
                public ChannelPipeline getPipeline() throws Exception {
                    return Channels.pipeline(new RemoteServerHandler());
                }
            });
            bootstrap.setOption("child.tcpNoDelay", true);
            bootstrap.setOption("child.keepAlive", true);
            bootstrap.setOption("child.connectTimeoutMillis", TimeUnit.SECONDS.toMillis(1));
            bootstrap.setOption("reuseAddress", true);
            bootstrap.setOption("tcpNoDelay", true);
            Channel channel = bootstrap.bind(new InetSocketAddress(hostname, port));
            return channel.getCloseFuture();

        }

        public RemoteActor register(Class< ? > actor, Object... constructorArgs) {
            return register(actor.getSimpleName(), actor, constructorArgs);
        }

        public RemoteActor register(String name, Class< ? > actor, Object... constructorArgs) {
            actorTable.put(name, actor);
            _logger.info("register an actor name:%s, class:%s", name, actor.getName());
            return new RemoteActor(name, actor, constructorArgs);
        }
    }

    public static class RemoteActor
    {
        private static Logger   _logger = Logger.getLogger(Actor.RemoteActor.class);

        public final String     name;
        public final Class< ? > actorClass;
        private Object          instance;

        public RemoteActor(String name, Class< ? > actorClass, Object... constructorArgs) {
            this.name = name;
            this.actorClass = actorClass;

            init(constructorArgs);
        }

        protected void init(Object... constructorArgs) {
            try {
                try {
                    Constructor< ? > constructor = actorClass.getConstructor(getParameterTypes(constructorArgs));
                    instance = constructor.newInstance(constructorArgs);
                }
                catch (NoSuchMethodException e) {
                    _logger.error(e);
                    instance = actorClass.newInstance();
                }
            }
            catch (Exception e) {
                _logger.error(e);
            }
        }

        public void call(String methodName, Object... args) {

            // TODO use client connection to send the function arguments to the remote server
            try {
                Method m = actorClass.getMethod(methodName, getParameterTypes(args));
                m.invoke(instance, args);
            }
            catch (Exception e) {
                _logger.error(e);
            }
        }

        Class< ? >[] getParameterTypes(Object... args) {
            Class< ? >[] argTypes = new Class< ? >[args.length];
            for (int i = 0; i < argTypes.length; ++i) {
                argTypes[i] = args[i].getClass();
            }
            return argTypes;
        }
    }

    private static class RemoteServerHandler extends SimpleChannelUpstreamHandler
    {

        @Override
        public void messageReceived(ChannelHandlerContext ctx, MessageEvent e) throws Exception {
            Channel ch = e.getChannel();
            // echo
            ch.write(e.getMessage());
        }

    }

    public static RemoteServer remote(String hostname, int port) {
        RemoteServer server = new RemoteServer(hostname, port);
        server.start();
        return server;
    }
}
