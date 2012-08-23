//--------------------------------------
//
// FASTAIndex.scala
// Since: 2012/06/22 5:09 PM
//
//--------------------------------------

package utgenome.weaver.lens

/**
 * Indexes of FASTA data
 *
 * @author leo
 */
trait FASTAIndex {

  def apply(seqName: String): DNASeq

  /**
   * Extract the sub sequence of the specified range [start, end).
   * @param seqName sequence name. e.g., chromosome name
   * @param start 0-based index (inclusive)
   * @param end 0-based index (exclusive)
   * @return
   */
  def subSequence(seqName: String, start: Int, end: Int): DNASeq

  def sequenceLength(seqName: String): Int

  def sequenceNames: Iterable[String]
}


class FASTAEntryIndex(val name: String, val description: String, val offset: Long, val length: Int)

trait FASTAIndexLike[Repr <: DNASeq with DNASeqOps[Repr]] {

  protected val seq : Repr
  protected val entry : Seq[FASTAEntryIndex]

  protected lazy val index: Map[String, FASTAEntryIndex] = (entry.map(e => e.name -> e)).toMap[String, FASTAEntryIndex]

  def sequenceLength(seqName: String) = index.apply(seqName).length

  def sequenceNames = index.keys

  /**
   * Retrieve a DNASeq of the given name
   * @param seqName
   * @return
   */
  def apply(seqName: String): DNASeq = new WrappedFASTASeq[Repr](seq, index(seqName).offset, sequenceLength(seqName))

  /**
   * Extract the sub sequence of the specified range [start, end).
   * @param seqName sequence name. e.g., chromosome name
   * @param start 0-based index (inclusive)
   * @param end 0-based index (exclusive)
   * @return
   */
  def subSequence(seqName: String, start: Int, end: Int): DNASeq = {
    val globalIndex = index(seqName).offset + start
    new WrappedFASTASeq[Repr](seq, globalIndex, end - start)
  }
}

class FASTAIndex2bit(protected val seq: ACGTSeq, protected val entry: Seq[FASTAEntryIndex])
  extends FASTAIndex
  with FASTAIndexLike[ACGTSeq] {

}

class FASTAIndex3bit(protected val seq: ACGTNSeq, protected val entry: Seq[FASTAEntryIndex])
  extends FASTAIndex
  with FASTAIndexLike[ACGTNSeq] {

}


/**
 * Wrapped DNASeq
 *
 * @param seq
 * @param offset
 * @param length
 */
class WrappedFASTASeq[Repr <: DNASeq with DNASeqOps[Repr]](seq: Repr, offset: Long, length: Long)
  extends DNASeq
  with DNASeqOps[Repr] {

  def domain = seq.domain

  def apply(index: Long) = seq(offset + index)

  def slice(start: Long, end: Long) = seq.slice(offset + start, offset + end)

  def numBases = length

  def count(base: DNA, start: Long, end: Long) = seq.count(base, offset + start, offset + end)

  /**
   * Count the number of occurrences of A, C, G and T letters within [start, end) at the same time.
   * This method is faster than repeating fastCount for each base.
   * @param start
   * @param end
   * @return
   */
  def count(start: Long, end: Long) = seq.count(offset + start, offset + end)

  /**
   * Create an materialized instance of this wrapped sequence
   * @return
   */
  protected def materialize = seq.slice(offset, offset+length)

  /**
   * Take the complement (not reversed) of this sequence
   * @return
   */
  def complement = materialize.complement

  /**
   * Create a reverse string of the this sequence. For example ACGT becomes TGCA
   *
   * @return Reverse sequence. The returned sequence is NOT a complement of
   *         the original sequence.
   */
  def reverse = materialize.reverse

  /**
   * Reverse complement of this sequence.
   * @return
   */
  def reverseComplement = materialize.reverseComplement
}