/*
 * Copyright 2012 Taro L. Saito
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, 
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.utgenome.weaver

import org.scalatest.FunSuite

//--------------------------------------
//
// LoggerTest.scala
// Since: 2012/01/04 13:28
//
//--------------------------------------

/**
 * Created at 2012/01/04 13:28
 * @author leo
 */
class LoggerTest extends FunSuite {

  class A extends Logger {

    def testSuite(f: LogFunction) {
       f.apply("Hello")
       f.apply("%s world! %d", "Hello", 2012)
    }

    val loggerType : List[LogFunction] = List(fatal, error, warn, info, debug, trace)
    loggerType.foreach { testSuite(_) }
  }

  test("display log messages in various log mode") {
    val a = new A
  }

}