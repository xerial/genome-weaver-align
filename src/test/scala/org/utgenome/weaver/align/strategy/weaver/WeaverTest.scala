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
package org.utgenome.weaver.align.strategy.weaver

import org.junit.Test
import scala.collection.JavaConversions._
import org.xerial.util.log.Logger
import org.scalatest.junit.{JUnitRunner, AssertionsForJUnit, ShouldMatchersForJUnit}
import org.scalatest.FunSuite
import org.scalatest.matchers.ShouldMatchers


class WeaverTest extends FunSuite with ShouldMatchers {

  val _logger: Logger = Logger.getLogger(classOf[WeaverTest])


  test("hello") {
    "hello world"
  }

  test("hello2") {

    val sb = new StringBuilder
    sb.append("ScalaTest is ")
    sb.append("fun!")
    sb.toString should be("ScalaTest is fun!")
    println(sb.toString)
  }

  sealed trait Fruit

  case class Apple() extends Fruit

  case class Banana() extends Fruit

  abstract class Control {
    def name: Unit
  }

  test("case class") {
    val e: Fruit = Apple()
    e match {
      case Apple() =>
      case Banana() =>
    }

    class A extends Control {
      def name {}
    }

    class B extends Control {
      def name {}
    }
  }


  test("get env") {
    for (e <- System.getenv().entrySet(); if e.getKey.toUpperCase.startsWith("TERM")) {
      _logger.debug("%s = %s".format(e.getKey, e.getValue.toString))
    }

  }

  test("pattern match") {

    class A(val name: String)
    class B(name: String) extends A(name)

    val a = new A("leo")
    val b = new B("yui")

    def fit(x: Any): String = {
      x match {
      case _:B => "this is B"
      case _:A => "this is A"
      case _ => "unknown"
    }
    }

    println(fit(a))
    println(fit(b))



  }

}

