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
trait FASTAIndex[A, +Repr] {

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

class FASTAIndex2bit(seq:ACGTSeq, entry:Seq[FASTAEntryIndex]) extends FASTAIndex[DNA2bit, DNASeq[DNA2bit]] {

  private lazy val index = (entry.map(e => e.name -> e)).toMap[String, FASTAEntryIndex]

  def apply(seqName: String) : DNASeq[DNA2bit] = new WrappedFASTASeq[DNA2bit, ACGTSeq](seq, index(seqName).offset, sequenceLength(seqName))

  /**
   * Extract the sub sequence of the specified range [start, end).
   * @param seqName sequence name. e.g., chromosome name
   * @param start 0-based index (inclusive)
   * @param end 0-based index (exclusive)
   * @return
   */
  def subSequence(seqName: String, start: Int, end: Int) : DNASeq[DNA2bit] = {
    val globalIndex = index(seqName).offset + start
    new WrappedFASTASeq[DNA2bit, ACGTSeq](seq, globalIndex, end - start)
  }

  def sequenceLength(seqName: String) = index.apply(seqName).length

  def sequenceNames = index.keys
}


/**
 * @param seq
 * @param offset
 * @param length
 */
class WrappedFASTASeq[+A, +Repr <: DNASeq[A] with DNASeqOps[A, Repr]](seq:Repr, offset:Long, length:Long)
  extends DNASeq[A]
  with DNASeqOps[A, Repr] {

  def domain = seq.domain

  def apply(index: Long) = seq(offset + index)

  def slice(start: Long, end: Long) = seq.slice(offset+start, offset+end)

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

  /**
   * Take the complement (not reversed) of this sequence
   * @return
   */
  def complement = seq.slice(offset, offset+length).complement

  /**
   * Create a reverse string of the this sequence. For example ACGT becomes TGCA
   *
   * @return Reverse sequence. The returned sequence is NOT a complement of
   *         the original sequence.
   */
  def reverse = seq.slice(offset, offset+length).reverse

  /**
   * Reverse complement of this sequence.
   * @return
   */
  def reverseComplement = seq.slice(offset, offset+length).reverseComplement
}