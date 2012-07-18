//--------------------------------------
//
// ACGTSeq.scala
// Since: 2012/06/19 3:42 PM
//
//--------------------------------------

package utgenome.weaver.lens

import java.util.Arrays

/**
 * Helper methods for managing 2bit encoding of DNA in an array of Long type
 * @author leo
 */
trait DNA2bit {

  def domain = DNA.exceptN

  // 2G (max of Java array size) * 8 (long byte size) * 8 (bit) / 2 (bit code) =  64G  (64G characters)
  val MAX_SIZE: Long = 2L * 1024L * 1024L * 1024L * 8L * 8L / 2L

  protected def minArraySize(numBases: Long): Int = {
    val bitSize: Long = numBases * 2L
    val blockBitSize: Long = 64L
    val arraySize: Long = (bitSize + blockBitSize - 1L) / blockBitSize
    if (arraySize > Integer.MAX_VALUE) {
      throw new IllegalArgumentException("Cannot create ACGTSequece more than %,d characters: %,d".format(MAX_SIZE, numBases))
    }
    arraySize.toInt
  }

  protected def blockIndex(basePos: Long): Int = (basePos >>> 5).toInt

  protected def blockOffset(basePos: Long): Int = (basePos & 0x1FL).toInt

}

/**
 * Utilities to build ACGTSeq
 */
object ACGTSeq {

  implicit def canBuildSeq : DNASeqBuilderFactory[ACGTSeq] = {
    new DNASeqBuilderFactory[ACGTSeq] {
      def apply() = new ACGTSeqBuilder()
    }
  }


  def newBuilder = new ACGTSeqBuilder()

  def newBuilder(sizeHint: Long) = new ACGTSeqBuilder()

  def apply(s: String): ACGTSeq = {
    val b = newBuilder
    b.sizeHint(s.length)
    for (ch <- s) {
      b += DNA(ch)
    }
    b.toACGTSeq
  }

  //  public static ACGTSeq loadFrom(File f) throws IOException {
  //    DataInputStream d = new DataInputStream(new BufferedInputStream(new FileInputStream(f), 4 * 1024 * 1024));
  //    try {
  //      return ACGTSeq.loadFrom(d);
  //    }
  //    finally {
  //      d.close();
  //    }
  //  }
  //
  //  public static ACGTSeq loadFrom(DataInputStream in) throws IOException {
  //    // The num bases must be always 2
  //
  //    long numBases = in.readLong();
  //    int longArraySize = minArraySize(numBases);
  //    long[] seq = new long[longArraySize];
  //    SnappyInputStream sin = new SnappyInputStream(in);
  //    int readBytes = sin.read(seq);
  //    return new ACGTSeq(seq, numBases);
  //  }
  //
  //  public void saveTo(DataOutputStream out) throws IOException {
  //    out.writeLong(this.numBases);
  //    SnappyOutputStream sout = new SnappyOutputStream(out);
  //    int longArraySize = minArraySize(this.numBases);
  //    sout.write(seq, 0, longArraySize);
  //    sout.flush();
  //  }
  //
  //  public void saveTo(File file) throws IOException {
  //    DataOutputStream d = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file)));
  //    try {
  //      saveTo(d);
  //    }
  //    finally {
  //      d.close();
  //    }
  //  }

}


/**
 * 2-bit encoded DNA Sequence of A, C, G and T. The maximum size this class can hold is
 * 2G (max of Java array size) * 8 (long byte size) * 8 (bit) / 2 (bit encoding) =  64G  (64G characters)
 *
 * To generate an instance of ACGTSeq, use ACGTSeq.newBuilder or ACGTSeq.apply
 *
 * @param seq 2-bit repre
 * @param numBases
 */
class ACGTSeq(private val seq: Array[Long], val numBases: Long)
  extends DNASeq
  with DNASeqOps[ACGTSeq]
  with DNA2bit
  with CharSequence {

  private var hash: Int = 0

  /**
   * Return the DNA base at the specified index
   * @param index
   * @return
   */
  def apply(index: Long): DNA = {
    // |---   64bit (32 characters) ---|
    // |ACGT|ACGT|ACGT| ...       |ACGT|

    val pos = blockIndex(index)
    val offset = blockOffset(index)
    val shift = offset << 1

    val code = (seq(pos) >>> shift) & 0x03
    DNA(code.toInt)
  }

  private def numFilledBlocks = (numBases / 32L).toInt

  private def lastBlock = seq.last & (~0L << ((numBases & 0x01FL).toInt * 2))

  override def hashCode() = {
    if (hash == 0) {
      var h = length * 31L
      var pos = 0
      for (i <- (0 until numFilledBlocks)) {
        h += seq(pos) * 31L
        pos += 1
      }
      val offset = (numBases % 32L).toInt
      if (offset > 0) {
        h += lastBlock * 31L
      }
      hash = h.toInt
    }
    hash
  }

  override def equals(obj: Any): Boolean = {
    obj match {
      case other: ACGTSeq => {
        if (this.numBases != other.numBases)
          false
        else {
          (0 until numFilledBlocks).find(i => this.seq(i) != other.seq(i)) match {
            case Some(x) => false
            case None => this.lastBlock == other.lastBlock
          }
        }
      }
      case _ => false
    }
  }

  /**
   * Create a reverse string of the this sequence. For example ACGT becomes TGCA
   *
   * @return Reverse sequence. The returned sequence is NOT a complement of
   *         the original sequence.
   */
  def reverse : ACGTSeq = {
    val b = new ACGTSeqBuilder(numBases)
    for(i <- 0L until numBases)
      b += apply(numBases - i - 1)
    b.result
  }


  /**
   * Take the complement (not reversed) of this sequence
   * @return
   */
  def complement: ACGTSeq = {
    val c = for (b <- seq) yield ~b
    new ACGTSeq(c, numBases)
  }

  def reverseComplement : ACGTSeq = {
    val b = new ACGTSeqBuilder(numBases)
    for(i <- 0L until numBases)
      b += apply(numBases - i - 1)
    val r = b.rawArray
    for(i <- 0 until seq.length) {
      r(i) = ~r(i)
    }
    new ACGTSeq(r, numBases)
  }


  protected def fastCount(v: Long, base: DNA): Long = {
    var r = ~0L
    r &= (if ((base.code & 0x02) == 0) ~v else v) >>> 1
    r &= (if ((base.code & 0x01) == 0) ~v else v)
    r &= 0x5555555555555555L
    // JVM optimizer is smart enough to replace this code to a pop count operation available in the CPU
    java.lang.Long.bitCount(r)
  }

  /**
   * Count the number of occurrences of the ACGT within the specified range [start, end)
   * @param base
   * @param start
   * @param end
   * @return the number of occurrences of the specified DNA base
   */
  def count(base: DNA, start: Long, end: Long): Long = {

    def countACGT: Long = {
      val sPos = blockIndex(start)
      val sOffset = blockOffset(start)

      val ePos = blockIndex(end - 1L)
      val eOffset = blockOffset(end - 1L)

      var count = 0L
      var numAsInMaskedRegion = 0L
      var pos = sPos
      while (pos <= ePos) {
        var mask: Long = ~0L
        if (pos == sPos) {
          mask <<= (sOffset << 1L)
          numAsInMaskedRegion += sOffset
        }
        if (pos == ePos) {
          val rMask = ~0L >>> (62L - (eOffset << 1))
          mask &= rMask
          numAsInMaskedRegion += 31L - eOffset
        }
        // Applying bit mask changes all bases in the masked region to As (code=00)
        val v: Long = seq(pos) & mask
        val popCount = fastCount(v, base)
        count += popCount
        pos += 1

      }
      if (base == DNA.A)
        count - numAsInMaskedRegion
      else
        count
    }

    base match {
      case DNA.N => 0L
      case _ => countACGT
    }

  }

  /**
   * Count the number of occurrences of A, C, G and T letters within [start, end) at the same time.
   * This method is faster than repeating fastCount for each base.
   * @param start
   * @param end
   * @return
   */
  def count(start: Long, end: Long): Array[Long] = {
    val count = Array.fill[Long](4)(0L)

    val sPos = blockIndex(start)
    val sOffset = blockOffset(start)

    val ePos = blockIndex(end - 1)
    val eOffset = blockOffset(end - 1)

    var numAsInMaskedRegion = 0L
    var pos = sPos
    while (pos <= ePos) {
      var mask: Long = ~0L
      if (pos == sPos) {
        mask <<= sOffset << 1L
        numAsInMaskedRegion += sOffset
      }
      if (pos == ePos) {
        val rMask = ~0L >>> (62L - (eOffset << 1))
        mask &= rMask
        numAsInMaskedRegion += 31L - eOffset
      }
      // Applying bit mask changes all bases in the masked region to As (code=00)
      val v: Long = seq(pos) & mask
      for (base <- DNA.exceptN) {
        count(base.code) += fastCount(v, base)
      }
      pos += 1
    }
    count(DNA.A.code) -= numAsInMaskedRegion
    count
  }

  /**
   * Extract a slice of the sequence [start, end)
   * @param start
   * @param end
   * @return
   */
  def slice(start: Long, end: Long): ACGTSeq = {
    if (start > end)
      sys.error("illegal argument start:%,d > end:%,d".format(start, end))

    val sliceLen = end - start
    val newSeq = new Array[Long](minArraySize(sliceLen))

    var i = 0L
    while (i < sliceLen) {
      val sPos = blockIndex(start + i)
      val sOffset = blockOffset(start + i)

      val dPos = blockIndex(i)
      val dOffset = blockOffset(i)

      var copyLen = 0L
      var l = 0L
      val v = seq(sPos) & (~0L << (sOffset * 2))
      if (sOffset == dOffset) {
        // no shift
        copyLen = 32L
        l = v
      }
      else if (sOffset < dOffset) {
        // left shift
        val shiftLen = dOffset - sOffset
        copyLen = 32L - dOffset
        l = v << (shiftLen * 2L)
      }
      else {
        // right shift
        val shiftLen = sOffset - dOffset
        copyLen = 32L - sOffset
        l = v >>> (shiftLen * 2L)
      }
      newSeq(dPos) |= l
      i += copyLen
    }

    new ACGTSeq(newSeq, sliceLen)
  }


}


/**
 * ACGT sequence builder
 * @param capacity the hint of number of bases to store
 */
class ACGTSeqBuilder(private var capacity: Long)
  extends DNA2bit
  with DNASeqBuilder[ACGTSeq]
{

  private var seq = new Array[Long](minArraySize(capacity))
  private var _numBases: Long = 0L

  /**
   * Create an empty sequence
   */
  def this() = this(32L)

  def numBases = _numBases

  def sizeHint(numBasesToStore: Long) {
    val arraySize = minArraySize(numBasesToStore)
    val newSeq = Arrays.copyOf(seq, arraySize)
    seq = newSeq
    capacity = arraySize * 32L
  }

  /**
   * Append a DNA base. Sequence other than A, C, G and T will be replaced to A
   * @param base
   */
  def +=(base: DNA): Unit = {
    val index = _numBases
    _numBases += 1
    if (index >= capacity) {
      val newCapacity = (index * 3L / 2L) + 64L
      sizeHint(newCapacity)
    }
    update(index, base)
  }


  /**
   * Set a DNA base at the specific index position. N will be replaced with A
   * @param index
   * @param base
   */
  def update(index: Long, base: DNA): Unit = {
    val pos = blockIndex(index)
    val offset = blockOffset(index)
    val shift = offset * 2

    // reset the target base. 3bit code N(code:100) will be trimmed to A (00)
    seq(pos) &= ~(0x3L << shift)
    seq(pos) |= (base.code & 0x03L) << shift
  }

  lazy val rawArray = {
    val size = minArraySize(_numBases)
    if (seq.length == size) seq else Arrays.copyOf(seq, size)
  }

  def result : ACGTSeq = toACGTSeq

  def toACGTSeq: ACGTSeq = new ACGTSeq(rawArray, _numBases)

  override def toString = result.toString

}
