//--------------------------------------
//
// OverlapSweeperTest.scala
// Since: 2012/07/03 2:59 PM
//
//--------------------------------------

package utgenome.weaver.lens

import xerial.silk.util.SilkSpec

/**
 * @author leo
 */
class OverlapSweeperTest extends SilkSpec {

  "OverlapSweeper" should {
    "report all overlapped intervals" in {

      // 1  2      5  6 7  8      10   12 13 14 15 16
      // |---------|              |----|  |-----|
      //    |-----------|         |----|      |-----|
      //              |---|

      val in = List(new Interval(1, 5), new Interval(2, 7), new Interval(6,8), new Interval(10, 12), new Interval(10, 12), new Interval(13, 15), new Interval(14, 16)).sorted(IntervalOrdering)
      val r = new OverlapSweeper(in)
      val overlapped = r.toSeq

      for(s <- overlapped) {
        for(c <- s.combinations(2)) {
          val a = c(0)
          val b = c(1)
          info("{%s, %s}", a, b)
          a.intersectWith(b) should be (true)
        }
      }

    }

  }

}