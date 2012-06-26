//--------------------------------------
//
// FASTATest.scala
// Since: 2012/06/22 10:19 AM
//
//--------------------------------------

package utgenome.weaver.lens

import xerial.silk.util.SilkSpec
import xerial.silk.util.io.TextDataProducer
import java.io.{StringReader, PrintWriter}
import util.Random
import org.xerial.util.FileResource
import org.scalatest.Tag


object ParserTest extends Tag("parser")
object ParserTest2 extends Tag("parser2")

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

    def randomFASTAReader = new StringReader(randomFASTA.readAsString)

    "parse fasta files" taggedAs(ParserTest) in {
      FASTA.read(randomFASTAReader) {
        stream =>
          for(r : FASTAEntryReader <- stream) {
            debug(r.name)
            for(line <- r.lines) {
              debug(line)
            }
          }
      }
    }

    def tgzFasta = FileResource.openByteStream(this.getClass, "sample-archive.fa.tar.gz")


    "parse tar.gz files"  taggedAs(ParserTest) in {
      tgzFasta should not be (null)
      FASTA.readTarGZ(tgzFasta) { stream =>
        for(r : FASTAEntryReader <- stream) {
          debug(r.name)
          for(line <- r.lines)
            debug(line)
        }
      }
    }


    "parse only the last chr"  taggedAs(ParserTest) in {
      FASTA.readTarGZ(tgzFasta) { stream =>
        for(r : FASTAEntryReader <- stream; if r.name == "chr3") {
          debug(r.name)
          debug(r.sequence)
        }
      }
    }


    "create index" taggedAs(ParserTest2) in {
      val index =  FASTA.create2bitIndexFromTarGZ(tgzFasta)
      val chrSet = index.sequenceNames.toSet
      List("chr1", "chr2", "chr3").forall(chrSet.contains(_)) should be (true)

      index("chr1")
    }


  }

  def randomFASTA() =
    new TextDataProducer {
      def produce(out: PrintWriter) {
        val r = new Random(0)
        for (i <- 0 until 3) {
          out.println(">chr%d".format(r.nextInt(21) + 1))
          for (w <- (0 until r.nextInt(1000)).sliding(80, 80)) {
            w.foreach(wi => out.print("ACGTN".charAt(r.nextInt(5))))
            out.println
          }
        }
      }
    }


}