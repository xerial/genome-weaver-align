//--------------------------------------
//
// FASTAIndex.scala
// Since: 2012/06/22 5:09 PM
//
//--------------------------------------

package xerial.silk.glens

/**
 * Indexes of FASTA data
 *
 * @author leo
 */
trait FASTAIndex[Repr <: DNASeq[Repr]] {

  def apply(seqName:String) : Repr

  /**
   * Extract the sub sequence of the specified range [start, end).
   * @param seqName sequence name. e.g., chromosome name
   * @param start 0-based index (inclusive)
   * @param end 0-based index (exclusive)
   * @return
   */
  def subSequence(seqName:String, start:Int, end:Int) : Repr

  def sequenceLength(seqName:String) : Int

  def sequenceNames : Iterable[String]
}


class FASTAEntryIndex(val name:String, val description:String, val offset:Long, val length:Int)

class FASTAIndex2bit(seq:ACGTSeq, entry:Seq[FASTAEntryIndex]) extends FASTAIndex[WrappedFASTASeq] {

  private lazy val index = (entry.map(e => e.name -> e)).toMap[String, FASTAEntryIndex]

  def apply(seqName: String) : WrappedFASTASeq = new WrappedFASTASeq(seq, index(seqName).offset, sequenceLength(seqName))

  /**
   * Extract the sub sequence of the specified range [start, end).
   * @param seqName sequence name. e.g., chromosome name
   * @param start 0-based index (inclusive)
   * @param end 0-based index (exclusive)
   * @return
   */
  def subSequence(seqName: String, start: Int, end: Int) : WrappedFASTASeq = {
    val globalIndex = index(seqName).offset + start
    new WrappedFASTASeq(seq, globalIndex, end - start)
  }

  def sequenceLength(seqName: String) = index.apply(seqName).length

  def sequenceNames = index.keys
}


/**
 * @param seq
 * @param offset
 * @param length
 */
class WrappedFASTASeq(seq:ACGTSeq, offset:Long, length:Long) extends DNASeq[WrappedFASTASeq] {

  def domain = seq.domain

  def apply(index: Long) = seq(offset + index)

  def slice(start: Long, end: Long) = new WrappedFASTASeq(seq, offset+start, end-start)

  def numBases = length

  def count(base: DNA, start: Long, end: Long) = seq.count(base, offset+start, offset+end)

  /**
   * Count the number of occurrences of A, C, G and T letters within [start, end) at the same time.
   * This method is faster than repeating fastCount for each base.
   * @param start
   * @param end
   * @return
   */
  def count(start: Long, end: Long) = seq.count(offset+start, offset+end)
}