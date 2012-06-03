//--------------------------------------
//
// ACGTSequenceTest.scala
// Since: 2012/06/03 10:50 AM
//
//--------------------------------------

package xerial.silk.glens

import xerial.silk.util.SilkSpec
import util.Random

/**
 * @author leo
 */
class ACGTSequenceTest extends SilkSpec {

  def compare(orig:String, seq:ACGTSequence) {
    orig.length should be (seq.length)
    seq.toACGTString should be (orig)
  }

  def randomSeq(len:Int) : String = {
    val r = new Random(len)
    val b = new StringBuilder
    for(i <- 0 until len)
      b += "ACGT".charAt(r.nextInt(4))
    b.result
  }

  "ACGTSequence" should {

    "construct instances from String" in {
      val seq = "AAACCGGTT"
      val s = ACGTSequence(seq)
      compare(seq, s)
    }

    "construct instances more than 32bp" in {
      val seq = randomSeq(1243)
      val w = ACGTSequence(seq)
      compare(seq, w)
    }

    "compute the reverse" in {
      val s = randomSeq(452)
      val c = ACGTSequence(s)
      c.numBases should be (s.length)
      val r = c.reverse

      compare(s.reverse, r)
    }

    "compute the reverse complement" in {
      def test(len:Int) {
        val s = randomSeq(len)
        val s_rc = {
          val sb = new StringBuilder
          for(i <- 0 until s.size reverse) {
            sb += DNA(s(i)).complement.toChar
          }
          sb.result
        }

        val c = ACGTSequence(s)
        val rc = c.reverseComplement
        compare(s_rc, rc)
      }

      test(3423)
      test(24)
    }

    "count base occurrences" in {

      def check(len:Int) {
        val s = randomSeq(len)
        val w = ACGTSequence(s)
        for(base <- DNA.exceptN) {
          for(x <- 0 until s.length; y <- x until s.length) {
            //debug("code:%s [%d, %d)", base, x, y)
            w.fastCount(base, x, y) should be (s.substring(x, y).count(c => c == base.toChar))
          }
        }
      }

      check(24)
      check(150)
    }

    "count base occurrences of ACGT at the same time" in {

      def check(len:Int) {
        val s = randomSeq(len)
        val w = ACGTSequence(s)
          for(x <- 0 until s.length; y <- x until s.length) {
            //debug("code:%s [%d, %d)", base, x, y)
            val count = w.fastCountACGT(x, y)
            for(base <- DNA.exceptN) {
              count(base.code) should be (s.substring(x, y).count(c => c == base.toChar))
            }
          }
      }
      check(24)
      check(150)
    }

  }

  "ACGTSequenceBuilder" should {
    "be capable to generate long DNA sequences" in {
      val seq = randomSeq(1000000)
      debug("generate ACGT sequence ...")
      val b = ACGTSequence.newBuilder
      for(slice <- seq.sliding(100, 100)) {
        b += slice
      }
      debug("done")
      val s = b.result
      compare(seq, s)
    }
  }

}