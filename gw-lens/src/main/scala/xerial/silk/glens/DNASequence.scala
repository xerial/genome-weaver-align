package xerial.silk.glens

import collection.mutable.ArrayBuffer
import xerial.silk.glens.DNA.N
import xerial.silk.util.{LogWriter, Logger}
import java.util.Arrays

/*
 * Copyright 2012 Taro L. Saito
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

//--------------------------------------
//
// DNASequence.scala
// Since: 2012/03/16 11:19
//
//--------------------------------------

/**
 * A base trait for classes representing DNA sequences
 *
 * @author leo
 */
trait DNASequence {

}

/**
 * Helper methods for managing 2bit encoding of DNA in an array of Long type
 * @author leo
 */
object DNA2bitEncoding {
  // 2G (max of Java array size) * 8 (long byte size) * 8 (bit) / 2 (bit code) =  64G  (64G characters)
  val MAX_SIZE: Long = 2L * 1024L * 1024L * 1024L * 8L * 8L / 2L

  def minArraySize(numBases: Long): Int = {
    val bitSize: Long = numBases * 2L
    val blockBitSize: Long = 64L
    val arraySize: Long = (bitSize + blockBitSize - 1L) / blockBitSize
    if (arraySize > Integer.MAX_VALUE) {
      throw new IllegalArgumentException("Cannot create ACGTSequece more than %,d characters: %,d".format(MAX_SIZE, numBases))
    }
    arraySize.toInt
  }

  def blockIndex(basePos: Long): Int = (basePos >>> 5).toInt


  def blockOffset(basePos: Long): Long = (basePos & 0x1FL) // This value must be Long to enable 64-bit shift using this value

}

/**
 * Utilities to build ACGTSeq
 */
object ACGTSeq {

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
 * To generate an instance of ACGTSeq, use {@link ACGTSeq$#newBuilder} or {@link ACGTSeq$#apply}.
 *
 * @param seq 2-bit repre
 * @param numBases
 */
class ACGTSeq(private val seq: Array[Long], val numBases: Long)
  extends DNASequence
  with Logger
  with CharSequence {

  import DNA2bitEncoding._

  private var hash: Int = 0

  private def longToIntCheck {
    if (numBases >= Integer.MAX_VALUE)
      sys.error("this method cannot be used when the sequence is larger than 2GB")
  }

  /**
   * length of this sequence
   * @return
   */
  def length = {
    longToIntCheck
    numBases.toInt
  }

  /**
   * Return string representation of this sequence
   * @return
   */
  def toACGTString: String = {
    val s = new StringBuilder
    var i = 0L
    while (i < numBases) {
      s += apply(i).toChar
      i += 1
    }
    s.result
  }

  override def toString = toACGTString

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
    val shift = offset << 1L

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
   * Take the complement (not reversed) of this sequence
   * @return
   */
  def complement: ACGTSeq = {
    val c = for (b <- seq) yield ~b
    new ACGTSeq(c, numBases)
  }

  /**
   * Create a reverse string of the this sequence. For example ACGT becomes TGCA
   *
   * @return Reverse sequence. The returned sequence is NOT a complement of
   *         the original sequence.
   */
  def reverse: ACGTSeq = {
    val buffer = new ACGTSeqBuilder(numBases)
    var i = 0L
    while (i < numBases) {
      buffer += apply(numBases - i - 1)
      i += 1
    }
    buffer.toACGTSeq
  }

  /**
   * Reverse complement of this sequence.
   * @return
   */
  def reverseComplement: ACGTSeq = {
    reverse.complement
  }


  private def fastCount(v: Long, base: DNA): Long = {
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
  def fastCount(base: DNA, start: Long, end: Long): Long = {

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
        //debug("pos:%d, popCount:%d, count:%d, numATotal:%d", pos, popCount, count, numAsInMaskedRegion)
        count += popCount
        pos += 1

      }
      if (base == DNA.A)
        count - numAsInMaskedRegion
      else
        count
    }

    base match {
      case N => 0L
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
  def fastCountACGT(start: Long, end: Long): Array[Long] = {
    val count = ArrayBuffer.fill[Long](4)(0L)

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
    count.toArray
  }

  def charAt(index: Int) = {
    apply(index).toChar
  }

  def subSequence(start: Int, end: Int) = slice(start, end)

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
  extends DNASequence {

  import DNA2bitEncoding._

  private var seq = new Array[Long](DNA2bitEncoding.minArraySize(capacity))
  private var _numBases: Long = 0L

  /**
   * Create an empty sequence
   */
  def this() = this(32L)

  def numBases = _numBases

  protected[glens] def sizeHint(numBasesToStore: Long) {
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

  def +=(seq: String): Unit = {
    for (ch <- seq) {
      this.+=(DNA(ch))
    }
  }

  /**
   * Set a DNA base at the specific index position. N will be replaced with A
   * @param index
   * @param base
   */
  def update(index: Long, base: DNA): Unit = {
    val pos = blockIndex(index)
    val offset = blockOffset(index)
    // Important: the shift length must be a Long value to enable 64bit-shift operations
    val shift: Long = offset * 2L

    // reset the target base. 3bit code N(code:100) will be trimmed to A (00)
    seq(pos) &= ~(0x3L << shift)
    seq(pos) |= (base.code & 0x03) << shift
  }

  def result = toACGTSeq

  def toACGTSeq: ACGTSeq = {
    val size = minArraySize(_numBases)
    val arr = if (seq.length == size) seq
    else {
      val newArr = Arrays.copyOf(seq, size)
      newArr
    }
    new ACGTSeq(arr, _numBases)
  }

  override def toString = result.toString

}