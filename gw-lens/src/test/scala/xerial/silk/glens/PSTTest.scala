package xerial.silk.glens

import xerial.silk.util.SilkSpec

//--------------------------------------
//
// PSTTest.scala
// Since: 2012/07/11 0:14
//
//--------------------------------------

/**
 * @author leo
 */
class PSTTest extends SilkSpec {

  "PST" should {
    "insert intervals" in {
      var p = new PST[Interval](20)
      val l = List(Interval(1, 5), Interval(2, 7),
        Interval(6,8), Interval(10, 12),
        Interval(10, 12), Interval(13, 15),
        Interval(14, 16))
      l.foreach(a => p = p.insert(a, a.start, a.end))
      p.size should be (l.length)
      debug(p)
    }
  }

}