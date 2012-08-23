package org.utgenome.weaver.align.strategy.weaver

import org.junit.runner.RunWith
import org.scalatest.junit.JUnitRunner
import org.xerial.util.log.Logger
import org.scalatest.matchers.ShouldMatchers
import collection.mutable.Stack
import org.scalatest.prop.PropertyChecks
import org.scalatest.{FlatSpec, Spec, BeforeAndAfter, FunSuite}


trait TestBlock {
  this: FunSuite =>
  def set {
    test("test1") {
      "hello"
    }
    test("test2") {

    }
  }

  def set2 {
    test("test3") {

    }
  }
}

@RunWith(classOf[JUnitRunner])
class SpecTest extends FunSuite with BeforeAndAfter with TestBlock with ShouldMatchers {

  val _logger = Logger.getLogger(classOf[SpecTest])


  def wrap(proc: => Any) {
    try {
      _logger.debug("begin")
      proc
    }
    finally {
      _logger.debug("end")
    }
  }


  before {
    _logger.debug("before test")
  }
  after {
    _logger.debug("after test")
  }

  test("addition") {
    val sum = 1 + 1
    assert(sum === 2)
  }

  test("subtraction") {
    val diff = 4 - 1
    assert(diff === 3)
    //reload
  }

  test("wrapping") {
    wrap {
      _logger.debug("hello world!")
    }

  }
  testsFor(set)
  testsFor(set2)


  test("should test") {
    val a = List(0, 1)
    a should have(size(2))
    a(0) should be(0)
    a(1) should be(1)
  }

  /*
  class PropTest extends FunSuite with PropertyChecks {


  }
  */

  def matchCode(x: Int) = {
    x match {
      case 1 => "Apple"
      case 2 => "Banana"
      case _ => "NonFruit"
    }
  }

  test("matcher matcher") {
    val s = matchCode(2)
    println(s)
  }

}


@RunWith(classOf[JUnitRunner])
class StackSpec extends Spec {

  describe("A Stack") {

    it("should pop values in last-in-first-out order") {
      val stack = new Stack[Int]
      stack.push(1)
      stack.push(2)
      assert(stack.pop() === 2)
      assert(stack.pop() === 1)
    }

    it("should throw NoSuchElementException if an empty stack is popped") {
      val emptyStack = new Stack[String]
      intercept[NoSuchElementException] {
        emptyStack.pop()
      }
    }
  }
}

class SampleSpec extends FlatSpec with ShouldMatchers {
  "A Stack name" should "pop values in last-in-first-out order" in {
    val stack = new Stack[Int]
    stack.push(1)
    stack.push(2)
    stack.pop() should equal (2)
    stack.pop() should equal (1)
  }
  it should "throw NoSuchElementException if an empty stack is popped" in {
    val emptyStack = new Stack[String]
    evaluating { emptyStack.pop() } should produce [NoSuchElementException]
  }
}
