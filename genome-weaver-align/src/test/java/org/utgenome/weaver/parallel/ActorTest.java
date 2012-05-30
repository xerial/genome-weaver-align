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
// ActorTest.java
// Since: 2011/10/19
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.parallel;

import org.junit.Test;
import org.utgenome.weaver.parallel.Actor.RemoteActor;
import org.utgenome.weaver.parallel.Actor.RemoteServer;
import org.xerial.util.log.Logger;

public class ActorTest
{
    private static Logger _logger = Logger.getLogger(ActorTest.class);

    public static class HelloActor
    {
        public void hello() {
            _logger.debug("Hello World!");
        }
    }

    @Test
    public void sample() throws Exception {
        RemoteServer remote = Actor.remote("localhost", 8990);
        RemoteActor remoteActor = remote.register(HelloActor.class);

        remoteActor.call("hello");
    }
}
