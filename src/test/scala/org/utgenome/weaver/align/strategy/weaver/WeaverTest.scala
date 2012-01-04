package org.utgenome.weaver.align.strategy.weaver
import org.junit.Test
import scala.collection.JavaConversions._
import org.xerial.util.log.Logger
import org.scalatest.junit.{JUnitRunner, AssertionsForJUnit, ShouldMatchersForJUnit}


class WeaverTest extends ShouldMatchersForJUnit {

  val _logger : Logger = Logger.getLogger(classOf[WeaverTest])

  @Test
  def hello {
    "hello world"
  }

  @Test
  def hello2 {

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

  @Test
  def caseClass {
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

  @Test
  def env {
    for(e <- System.getenv().entrySet(); if e.getKey.toUpperCase.startsWith("TERM")) {
       _logger.debug("%s = %s".format(e.getKey, e.getValue.toString))
    }
    
  }


}

