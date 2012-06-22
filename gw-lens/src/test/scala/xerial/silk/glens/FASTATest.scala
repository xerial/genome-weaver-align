//--------------------------------------
//
// FASTATest.scala
// Since: 2012/06/22 10:19 AM
//
//--------------------------------------

package xerial.silk.glens

import xerial.silk.util.SilkSpec

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
   }


}