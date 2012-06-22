//--------------------------------------
//
// FASTATest.scala
// Since: 2012/06/22 10:19 AM
//
//--------------------------------------

package xerial.silk.glens

import xerial.silk.util.SilkSpec
import xerial.silk.util.io.TextDataProducer
import java.io.PrintWriter
import util.Random

/**
 * @author leo
 */
class FASTATest extends SilkSpec {

  "FASTA" should {
    "extract sequence name" in {
      val c1 = FASTA.extractSequenceNameFrom(">chr10")
      c1 should be("chr10")

      val c2 = FASTA.extractSequenceNameFrom("> chrX human chromosome X")
      c2 should be("chrX")

      val c3 = FASTA.extractSequenceNameFrom(">seqname|separated by |hirozotal bar")
      c3 should be("seqname")

    }

    "create random seq" in {
      val f = randomFASTA
      debug(f)
    }

  }

  def randomFASTA() = {
    val p = new TextDataProducer {
      def produce(out: PrintWriter) {
        val r = new Random(0)
        for (i <- 0 until 3) {
          out.println(">chr%d".format(r.nextInt(21) + 1))
          for(w <- (0 until r.nextInt(1000)).sliding(80, 80)) {
            w.foreach(wi => out.print("ACGTN".charAt(r.nextInt(5))))
            out.println
          }
        }
      }
    }
    p.lines.mkString("\n")
  }


}