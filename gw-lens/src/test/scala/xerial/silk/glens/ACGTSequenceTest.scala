//--------------------------------------
//
// ACGTSequenceTest.scala
// Since: 2012/06/03 10:50 AM
//
//--------------------------------------

package xerial.silk.glens

import xerial.silk.util.SilkSpec

/**
 * @author leo
 */
class ACGTSequenceTest extends SilkSpec {

  "ACGTSequence" should {

    "construct instances from String" in {
      val s = ACGTSequence("AAACCGGTT")
      s should have length (9)



    }


  }


}