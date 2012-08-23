//--------------------------------------
//
// ACGTNSeq.scala
// Since: 2012/06/19 4:31 PM
//
//--------------------------------------

package utgenome.weaver.lens

import java.util
import util.Arrays


/**
 * Helper methods for managing 2bit encoding of DNA in an array of Long type
 * @author leo
 */
trait DNA3bit {

  def domain = DNA.values

  // 2G (max of Java array size) * 8 (long byte size) * 8 (bit) / 3 (bit code) =  42G  (64G characters)
  val MAX_SIZE: Long = 2L * 1024L * 1024L * 1024L * 8L * 8L / 3L

  protected def minArraySize(numBases: Long): Int = {
    val totalBitSize: Long = numBases * 3L
    val blockBitSize: Long = 64L * 3L // 3 long blocks are one set
    val arraySize: Long = ((totalBitSize + blockBitSize - 1L) / blockBitSize) * 3L
    if (arraySize > Integer.MAX_VALUE) {
      throw new IllegalArgumentException("Cannot create ACGTSequece more than %,d characters: %,d".format(MAX_SIZE, numBases))
    }
    arraySize.toInt
  }

  protected def blockIndex(basePos: Long): Int = (basePos >>> 6).toInt

  protected def blockOffset(basePos: Long): Int = (basePos & 0x03FL).toInt

  protected def numFilledBlocks(numBases: Long) = (numBases / 64L * 3L).toInt

  protected def lookupBase(seq:Array[Long], index:Long) : DNA = {
    val pos: Int = blockIndex(index)
    val offset: Int = blockOffset(index)
    val shift = 62 - ((index & 0x1FL) << 1).toInt

    val nFlag: Long = seq(pos * 3) & (1L << (63 - offset))
    val code: Int = (seq(pos * 3 + (offset >> 5) + 1) >>> shift).toInt & 0x03
    if (nFlag == 0) DNA(code) else DNA.N
  }

  protected def updateBase(seq:Array[Long], index:Long, base:DNA)  {
    val pos = blockIndex(index)
    val offset = blockOffset(index)
    val shift = (offset & 0x1F) << 1

    val code = base.code
    seq(pos * 3) &= ~(1L << (63 - offset))
    seq(pos * 3) |= ((code >>> 2) & 0x01L) << (63 - offset)
    val bPos: Int = pos * 3 + (offset >> 5) + 1
    seq(bPos) &= ~(0xC000000000000000L >>> shift)
    seq(bPos) |= (code & 0x03L) << (62 - shift)
  }

}


object ACGTNSeq {

  def newBuilder = new ACGTNSeqBuilder
  def newBuilder(numBases:Long) = new ACGTNSeqBuilder(numBases)

  def apply(s:String) = {
    val b = newBuilder(s.length)
    for(ch <- s) {
      b += DNA(ch)
    }
    b.result
  }

  implicit def canBuildSeq : DNASeqBuilderFactory[ACGTNSeq] = {
    new DNASeqBuilderFactory[ACGTNSeq] {
      def apply() = new ACGTNSeqBuilder()
    }
  }


}

/**
 * 3-bit encoding of DNA sequence, representing A, C, G, T and N.
 *
 * This sequence holds ACGTN using a sequence of three Long blocks.
 * The first block is for Ns  (0 or 1). The second and third blocks hold 2-bit encoding of ACGT.
 * |N0 ... N64|B0 ... B31|B32 ... B63|
 *
 *
 * @author leo
 */
class ACGTNSeq(private val seq: Array[Long], val numBases: Long)
  extends DNASeq
  with DNASeqOps[ACGTNSeq]
  with DNA3bit {


  protected var hash: Int = 0

  override def hashCode = {
    if (hash == 0) {
      var h = numBases * 31L
      var pos = 0
      val filledBlockSize = numFilledBlocks(numBases)
      while (pos < filledBlockSize) {
        h += seq(pos) * 31L
        pos += 1
      }
      val offset = (numBases % 64L).toInt
      if (offset > 0) {
        h += (seq(pos) & (~0L << (64 - offset))) * 31L
        h += (seq(pos + 1) & (if (offset < 32) ~0L << ((32 - offset) * 2) else ~0L)) * 31L
        h += (seq(pos + 2) & (if (offset < 32) 0L else ~0L << (64 - offset) * 2)) * 31L
      }
      hash = h.toInt
    }
    hash
  }

  override def equals(obj: Any) = {

    def eligible: Option[ACGTNSeq] = {
      if (!(obj.isInstanceOf[ACGTNSeq]))
        None
      else
        Some(obj.asInstanceOf[ACGTNSeq])
    }

    def hasSameFilledBlocks(other: ACGTNSeq) =
      (0 until numFilledBlocks(numBases)).forall(i => this.seq(i) == other.seq(i))

    val isEqual = eligible.filter(this.numBases == _.numBases).filter(hasSameFilledBlocks).filter {
      other =>
        val offset: Int = (numBases % 64L).toInt
        if (offset == 0)
          true
        else {
          val mask: Array[Long] = new Array[Long](3)
          mask(0) = ~0L << 64 - offset
          mask(1) = if (offset < 32) ~0L << (32 - offset) * 2 else ~0L
          mask(2) = if (offset <= 32) 0L else ~0L << (64 - offset) * 2
          val pos = numFilledBlocks(numBases)
          (0 until mask.length).forall(i => (seq(pos + i) & mask(i)) == (other.seq(pos + i) & mask(i)))
        }
    }
    isEqual.isDefined
  }


  def apply(index: Long) = lookupBase(seq, index)

  def slice(start: Long, end: Long) = {
    if (start > end)
      sys.error("invalid range [%d, %d)".format(start, end))

    val len = end - start
    val dest = new Array[Long](minArraySize(len))

    var i = 0L
    while (i < len) {
      val sPos = blockIndex(start + i)
      val sOffset = blockOffset(start + i)
      val dPos = blockIndex(i)
      val dOffset = blockOffset(i)

      var copyLen = 0L
      var n = seq(sPos * 3)
      var h = seq(sPos * 3 + 1)
      var l = seq(sPos * 3 + 2)
      if (sOffset == dOffset) {
        copyLen = 64L
      }
      else if (sOffset < dOffset) {
        // right shift
        val shiftLen = dOffset - sOffset;
        copyLen = 64 - dOffset
        // Copy Ns
        n >>>= shiftLen
        // Copy ACGT blocks
        if (shiftLen < 32) {
          l = (h << (64 - shiftLen * 2)) | (l >>> shiftLen * 2)
          h >>>= shiftLen * 2
        }
        else {
          l = h >>> (shiftLen - 32) * 2
          h = 0L
        }
      }
      else {
        // left shift
        val shiftLen = sOffset - dOffset
        copyLen = 64 - sOffset
        // Copy Ns
        n <<= shiftLen
        // Copy ACGT blocks
        if (shiftLen < 32) {
          h = (h << shiftLen * 2) | (l >>> (64 - shiftLen * 2))
          l <<= shiftLen * 2
        }
        else {
          h = l << (shiftLen - 32) * 2
          l = 0L
        }
      }
      dest(dPos * 3) |= n
      dest(dPos * 3 + 1) |= h
      dest(dPos * 3 + 2) |= l

      i += copyLen
    }

    new ACGTNSeq(dest, len)
  }

  def count(base: DNA, start: Long, end: Long) = {
    var count: Long = 0
    if (base == DNA.N) {
      var sPos: Int = (start >>> 6).toInt
      var sOffset: Int = (start & 0x3FL).toInt
      val ePos: Int = ((end + 64L - 1L) >>> 6).toInt
      while (sPos < ePos) {
        {
          var mask: Long = ~0L
          if (sOffset != 0) {
            mask >>>= sOffset
            sOffset = 0
          }
          if (sPos == ePos - 1) {
            val eOffset: Int = (end & 0x3FL).toInt
            val rMask: Long = if ((eOffset == 0)) ~0L else ~((1L << (64 - eOffset)) - 1)
            mask &= rMask
          }
          count += java.lang.Long.bitCount(seq(sPos * 3) & mask)
        }
        sPos += 1
      }
    }
    else {
      var sPos: Int = (start >>> 5).toInt
      var sOffset: Int = (start & 0x1FL).toInt
      val ePos: Int = ((end + 32L - 1L) >>> 5).toInt
      while (sPos < ePos) {
        {
          var mask: Long = ~0L
          if (sOffset != 0) {
            mask >>>= sOffset * 2
            sOffset = 0
          }
          val bIndex: Int = sPos / 2 * 3
          val block: Int = sPos % 2
          val v: Long = seq(bIndex + 1 + block)
          val nFlag: Long = interleave32With0(seq(bIndex) >>> (32 * (1 - block)))
          if (sPos == ePos - 1) {
            val eOffset: Int = (end & 0x1FL).toInt
            val rMask: Long = if ((eOffset == 0)) ~0L else ~((1L << (32 - eOffset) * 2) - 1)
            mask &= rMask
          }
          var r: Long = ~0L
          r &= (if ((base.code & 0x02) == 0) ~v else v) >>> 1
          r &= (if ((base.code & 0x01) == 0) ~v else v)
          r &= 0x5555555555555555L
          r &= ~nFlag
          r &= mask
          count += java.lang.Long.bitCount(r)
        }
        sPos += 1
      }
    }
    count
  }

  /**
   * Count the number of occurrences of A, C, G and T letters within [start, end) at the same time.
   * This method is faster than repeating fastCount for each base.
   * @param start
   * @param end
   * @return
   */
  def count(start: Long, end: Long) = {
    val count = new Array[Long](5)

    // Count A, C, G, T
    var sPos = (start >>> 5).toInt
    var sOffset = (start & 0x1FL).toInt
    val ePos = ((end + 32L - 1L) >>> 5)

    while (sPos < ePos) {
      var mask = ~0L
      if (sOffset != 0) {
        mask >>>= sOffset * 2
        sOffset = 0
      }
      val bIndex = sPos / 2 * 3
      val block = sPos % 2
      val v = seq(bIndex + 1 + block)
      val nFlag = interleave32With0(seq(bIndex) >>> (32 * (1 - block)))
      if (sPos == ePos - 1) {
        val eOffset = (end & 0x1FL).toInt
        val rMask = if (eOffset == 0) ~0L else ~((1L << (32 - eOffset) * 2) - 1)
        mask &= rMask
      }

      // count A, C, G and T
      for (base <- DNA.exceptN) {
        var r = ~0L
        r &= (if ((base.code & 0x02) == 0) ~v else v) >>> 1
        r &= (if ((base.code & 0x01) == 0) ~v else v)
        r &= 0x5555555555555555L
        r &= ~nFlag
        r &= mask
        count(base.code) += java.lang.Long.bitCount(r)
      }

      // count N
      count(DNA.N.code) += java.lang.Long.bitCount(nFlag & mask)

      sPos += 1
    }

    count
  }

  /**
   * Interleave low 32bits (in a long value) with 0s. For example, 11110011 (8
   * bit value) becomes 0101010100000101 (16 bit value)
   *
   * @param input
   * @return
   */
  private def interleave32With0(input: Long): Long = {
    var v = input
    v = ((v & 0xFFFF0000L) << 16) | (v & 0x0000FFFFL)
    v = ((v << 8) | v) & 0x00FF00FF00FF00FFL
    v = ((v << 4) | v) & 0x0F0F0F0F0F0F0F0FL
    v = ((v << 2) | v) & 0x3333333333333333L
    v = ((v << 1) | v) & 0x5555555555555555L
    v
  }

  private def complement(c: Array[Long]) : Array[Long] = {
    val numBlocks = seq.length / 3
    for (i <- 0 until numBlocks) {
      c(i * 3) = c(i * 3)
      c(i * 3 + 1) = ~(c(i * 3 + 1))
      c(i * 3 + 2) = ~(c(i * 3 + 2))
    }
    c
  }


  private def rawReverse = {
    val b = new ACGTNSeqBuilder(numBases)
    for (i <- (numBases - 1L) to 0L by -1L) {
      b += apply(i)
    }
    b.rawArray
  }

  def reverse = {
    new ACGTNSeq(rawReverse, this.numBases)
  }


  def complement : ACGTNSeq = {
    val c = new Array[Long](seq.length)
    new ACGTNSeq(complement(c), this.numBases)
  }


  def reverseComplement: ACGTNSeq = {
    val r = rawReverse
    new ACGTNSeq(complement(r), this.numBases)
  }
}


class ACGTNSeqBuilder(private var capacity:Long)
  extends DNASeqBuilder[ACGTNSeq]
  with DNA3bit
{

  def this() = this(10L)

  private var seq = new Array[Long](minArraySize(capacity))
  private var _numBases = 0L

  def numBases = _numBases

  capacity = minArraySize(capacity) / 3L * 64L


  def sizeHint(newSize: Long) {
    if (capacity <= newSize) {
      val requiredArraySize = minArraySize(newSize)
      if (requiredArraySize >= seq.length) {
        val newArraySize = minArraySize((newSize * 1.5 + 64L).toLong)
        val newArray = util.Arrays.copyOf(seq, newArraySize)
        capacity = newArraySize / 3L * 64L
        seq = newArray
      }
    }
  }

  def +=(base: DNA) {
    val index = _numBases
    _numBases += 1

    // |N0 ... N63|B0 B1 ....  B31|B32 B33 ... B63|
    sizeHint(index)
    updateBase(seq, index, base)
  }

  lazy val rawArray = {
    val size = minArraySize(_numBases)
    if (seq.length == size) seq else Arrays.copyOf(seq, size)
  }

  def result = new ACGTNSeq(rawArray, _numBases)

}


/**
 * Mutable version of ACGTNSeq
 * @param seq
 * @param numBases
 */
class ACGTNSeqBuffer(private val seq: Array[Long], override val numBases: Long)
  extends ACGTNSeq(seq, numBases) {

  def update(index: Long, base: DNA) {
    updateBase(seq, index, base)
  }

}