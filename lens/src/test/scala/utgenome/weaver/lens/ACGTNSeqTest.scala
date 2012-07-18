//--------------------------------------
//
// ACGTNSeqTest.scala
// Since: 2012/06/21 4:58 PM
//
//--------------------------------------

package utgenome.weaver.lens

import xerial.silk.util.SilkSpec
import util.Random

/**
 * @author leo
 */
class ACGTNSeqTest extends SilkSpec {

  def compare(orig: String, seq: ACGTNSeq) {
    orig.length should be(seq.length)
    seq.toACGTString should be(orig)
  }

  def randomSeq(len: Int): String = {
    val r = new Random(len)
    val b = new StringBuilder
    for (i <- 0 until len)
      b += "ACGTN".charAt(r.nextInt(5))
    b.result
  }

  "ACGTNSeq (3-bit encoding)" should {

    "construct instances from String" in {
      val seq = "AAACCGGTT"
      val s = ACGTNSeq(seq)
      compare(seq, s)
    }

    "construct instances more than 32bp" in {
      val seq = randomSeq(1243)
      val w = ACGTNSeq(seq)
      compare(seq, w)
    }

    "compute the reverse" in {
      val s = randomSeq(452)
      val c = ACGTNSeq(s)
      c.numBases should be(s.length)
      val r = c.reverse

      compare(s.reverse, r)
    }

    "compute the reverse complement" in {
      def test(len: Int) {
        val s = randomSeq(len)
        val s_rc = {
          val sb = new StringBuilder
          for (i <- 0 until s.size reverse) {
            sb += DNA(s(i)).complement.toChar
          }
          sb.result
        }

        val c = ACGTNSeq(s)
        val rc = c.reverseComplement
        compare(s_rc, rc)
      }

      test(3423)
      test(24)
    }

    "slice sub strings" in {
      def test(len: Int) {
        val s = randomSeq(len)
        val w = ACGTNSeq(s)
        for (x <- 0 until s.length; y <- x until s.length) {
          w.slice(x, y).toACGTString should be(s.substring(x, y))
        }
      }

      test(24)
      test(100)
    }

    "count base occurrences" in {

      def check(len: Int) {
        val s = randomSeq(len)
        val w = ACGTNSeq(s)
        for (base <- DNA.exceptN) {
          for (x <- 0 until s.length; y <- x until s.length) {
            //debug("code:%s [%d, %d)", base, x, y)
            w.count(base, x, y) should be(s.substring(x, y).count(c => c == base.toChar))
          }
        }
      }

      check(24)
      check(150)
    }

    "count base occurrences of ACGT at the same time" in {

      def check(len: Int) {
        val s = randomSeq(len)
        val w = ACGTNSeq(s)
        for (x <- 0 until s.length; y <- x until s.length) {
          //debug("code:%s [%d, %d)", base, x, y)
          val count = w.count(x, y)
          for (base <- DNA.exceptN) {
            count(base.code) should be(s.substring(x, y).count(c => c == base.toChar))
          }
        }
      }
      check(24)
      check(150)
    }

  }

  "ACGTNSeqBuilder (3-bit)" should {
    "be capable to generate long DNA sequences" in {
      val seq = randomSeq(1000000)
      val b = ACGTNSeq.newBuilder
      for (slice <- seq.sliding(100, 100)) {
        b += slice
      }
      val s = b.result
      val cmp = seq == s.toACGTString
      cmp should be (true)
    }
  }

}