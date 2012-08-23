/*--------------------------------------------------------------------------
 *  Copyright 2012 utgenome.org
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
//
// SpecTest2.scala
// Since: 2012/01/03 12:17
//
//--------------------------------------
package org.utgenome.weaver
package align.strategy.weaver

import scala.actors.Actor._
import org.scalatest.FunSuite

/**
 * Created at 2012/01/03 12:17
 * @author leo
 */
class SpecTest2 extends FunSuite {

  test("actor test") {
    actor {
      loopWhile(true) {

      }
    }

  }


}