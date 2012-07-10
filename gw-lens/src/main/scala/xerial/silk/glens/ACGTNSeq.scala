//--------------------------------------
//
// ACGTNSeq.scala
// Since: 2012/06/19 4:31 PM
//
//--------------------------------------

package xerial.silk.glens


/**
 * Helper methods for managing 2bit encoding of DNA in an array of Long type
 * @author leo
 */
trait DNA3bitEncoding {
  // 2G (max of Java array size) * 8 (long byte size) * 8 (bit) / 3 (bit code) =  42G  (64G characters)
  val MAX_SIZE: Long = 2L * 1024L * 1024L * 1024L * 8L * 8L / 3L

  def minArraySize(numBases: Long): Int = {
    val bitSize: Long = numBases * 3L
    val blockBitSize: Long = 64L * 3L // 3 long blocks are one set
    val arraySize: Long = (bitSize + blockBitSize - 1L) / blockBitSize
    if (arraySize > Integer.MAX_VALUE) {
      throw new IllegalArgumentException("Cannot create ACGTSequece more than %,d characters: %,d".format(MAX_SIZE, numBases))
    }
    arraySize.toInt
  }

  def blockIndex(basePos: Long): Int = (basePos >>> 6).toInt


  def blockOffset(basePos: Long): Long = (basePos & 0x3FL) // This value must be Long to enable 64-bit shift using this value

  def numFilledBlocks(numBases: Long) = (numBases / 64L * 3L).toInt

}


object ACGTNSeq {


}

/**
 * 3-bit encoding of DNA sequence, representing A, C, G, T and N
 *
 * @author leo
 */
class ACGTNSeq(private val seq: Array[Long], val numBases: Long)
  extends DNASeq[ACGTNSeq]
  with DNA3bitEncoding {

  lazy private val hash: Int = {
    var h = numBases * 31L
    var pos = 0
    val filledBlockSize = numFilledBlocks(numBases)
    while (pos < filledBlockSize) {
      h += seq(pos) * 31L
      pos += 1
    }
    val offset = (numBases % 64L).toInt
    if (offset > 0) {
      h += (seq(pos) & (~0L << 64 - offset)) * 31L
      h += (seq(pos + 1) & (if (offset < 32) ~0L << (32 - offset) * 2 else ~0L)) * 31L
      h += (seq(pos + 2) & (if (offset < 32) 0L else ~0L << (64 - offset) * 2)) * 31L
    }
    h.toInt
  }

  override def hashCode = hash

//  override def equals(obj: Any) = {
//
//    def cast = {
//      if (!(obj.isInstanceOf[ACGTNSeq]))
//        None
//      else
//        Some(obj.asInstanceOf[ACGTNSeq])
//    }
//
//    cast.filter(other => this.numBases == other.numBases)
//
//
//
//
//    val isEqual = eligible.flatMap {
//      other =>
//
//        val filledBlockSize: Int = numFilledBlocks(numBases)
//        var pos: Int = 0
//        while (pos < filledBlockSize) {
//          if (this.seq(pos) != other.seq(pos))
//            return false
//          pos += 1
//        }
//
//        val offset: Int = (numBases % 64L).toInt
//
//        if (offset > 0) {
//          val mask: Array[Long] = new Array[Long](3)
//          mask(0) = ~0L << 64 - offset
//          mask(1) = if (offset < 32) ~0L << (32 - offset) * 2 else ~0L
//          mask(2) = if (offset <= 32) 0L else ~0L << (64 - offset) * 2
//          var i: Int = 0
//          while (i < mask.length) {
//            if ((seq(pos + i) & mask(i)) != (other.seq(pos + i) & mask(i)))
//              return false
//            i += 1
//          }
//        }
//        true
//    }
//
//    isEqual.
//  }


  def apply(index: Long) = {
    val pos: Int = blockIndex(index)
    val offset: Int = blockOffset(index).toInt
    val shift = 62L - ((index & 0x1FL) << 1)

    val nFlag: Long = seq(pos * 3) & (1L << (63 - offset))
    val code: Int = (seq(pos * 3 + (offset >> 5) + 1) >>> shift).toInt & 0x03
    if (nFlag == 0) DNA(code) else DNA.N
  }

  def slice(start: Long, end: Long) = null

  def count(base: DNA, start: Long, end: Long) = 0L

  /**
   * Count the number of occurrences of A, C, G and T letters within [start, end) at the same time.
   * This method is faster than repeating fastCount for each base.
   * @param start
   * @param end
   * @return
   */
  def count(start: Long, end: Long) = null

}