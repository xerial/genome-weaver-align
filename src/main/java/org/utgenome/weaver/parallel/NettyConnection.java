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
// NettyConnection.java
// Since: 2011/10/18
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.parallel;

import java.net.InetSocketAddress;
import java.net.SocketAddress;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.jboss.netty.bootstrap.ClientBootstrap;
import org.jboss.netty.bootstrap.ServerBootstrap;
import org.jboss.netty.channel.Channel;
import org.jboss.netty.channel.ChannelFuture;
import org.jboss.netty.channel.ChannelHandlerContext;
import org.jboss.netty.channel.ChannelPipeline;
import org.jboss.netty.channel.ChannelPipelineFactory;
import org.jboss.netty.channel.ChannelStateEvent;
import org.jboss.netty.channel.Channels;
import org.jboss.netty.channel.ExceptionEvent;
import org.jboss.netty.channel.MessageEvent;
import org.jboss.netty.channel.SimpleChannelUpstreamHandler;
import org.jboss.netty.channel.StaticChannelPipeline;
import org.jboss.netty.channel.socket.nio.NioClientSocketChannelFactory;
import org.jboss.netty.channel.socket.nio.NioServerSocketChannelFactory;
import org.jboss.netty.handler.codec.serialization.ObjectDecoder;
import org.jboss.netty.handler.codec.serialization.ObjectEncoder;
import org.xerial.util.log.Logger;

/**
 * An interface to transfer messages through Netty
 * 
 * @author leo
 * 
 */
public class NettyConnection
{
    private static Logger _logger = Logger.getLogger(NettyConnection.class);

    private static class ActiveRemoteClientPipelineFactory implements ChannelPipelineFactory
    {

        @Override
        public ChannelPipeline getPipeline() throws Exception {

            return new StaticChannelPipeline(new ObjectEncoder(), new ObjectDecoder(), new ActiveRemoteClientHandler());
        }

    }

    private static class ActiveRemoteClientHandler extends SimpleChannelUpstreamHandler
    {
        @Override
        public void messageReceived(ChannelHandlerContext ctx, MessageEvent e) throws Exception {
            _logger.info("Client: recieved %s", e.getMessage());
        }
    }

    public static boolean connect(String hostname, int port) {
        return connect(new InetSocketAddress(hostname, port));
    }

    public static boolean connect(SocketAddress remoteAddress) {

        ClientBootstrap bootstrap = new ClientBootstrap(new NioClientSocketChannelFactory(
                Executors.newCachedThreadPool(), Executors.newCachedThreadPool()));
        bootstrap.setPipelineFactory(new ActiveRemoteClientPipelineFactory());
        bootstrap.setOption("tcpNoDelay", true);
        bootstrap.setOption("keepAlive", true);

        ChannelFuture connection = bootstrap.connect(remoteAddress);
        Channel openChannel = connection.awaitUninterruptibly().getChannel();

        if (!connection.isSuccess()) {
            _logger.error("connection failed (remote address:%s): %s", remoteAddress, connection.getCause());
            return false;
        }

        return true;
    }

    private static class RemoteServer extends SimpleChannelUpstreamHandler
    {
        @Override
        public void channelOpen(ChannelHandlerContext ctx, ChannelStateEvent e) throws Exception {
            _logger.info("Server: channel opened");
        }

        @Override
        public void messageReceived(ChannelHandlerContext ctx, MessageEvent e) throws Exception {
            _logger.info("Server: recieved %s", e.getMessage());

        }

        @Override
        public void exceptionCaught(ChannelHandlerContext ctx, ExceptionEvent e) throws Exception {
            _logger.info("error: %s", e.getCause());
        }
    }

    private static class RemoteServerPipelineFactory implements ChannelPipelineFactory
    {
        @Override
        public ChannelPipeline getPipeline() throws Exception {
            return Channels.pipeline(new ObjectEncoder(), new ObjectDecoder(), new RemoteServer());
        }

    }

    public static ChannelFuture launchServer(String hostname, int port) {
        return launchServer(new InetSocketAddress(hostname, port));
    }

    public static ChannelFuture launchServer(SocketAddress serverAddress) {
        ServerBootstrap bootstrap = new ServerBootstrap(new NioServerSocketChannelFactory(
                Executors.newCachedThreadPool(), Executors.newCachedThreadPool()));

        bootstrap.setPipelineFactory(new RemoteServerPipelineFactory());
        bootstrap.setOption("child.tcpNoDelay", true);
        bootstrap.setOption("child.keepAlive", true);
        bootstrap.setOption("child.connectTimeoutMillis", TimeUnit.SECONDS.toMillis(1));
        bootstrap.setOption("reuseAddress", true);
        bootstrap.setOption("tcpNoDelay", true);
        Channel channel = bootstrap.bind(serverAddress);
        return channel.getCloseFuture();
    }
}
