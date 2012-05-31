package xerial.silk.glens

import collection.mutable.ArrayBuffer
import xerial.silk.glens.DNA.N

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
 * @author leo
 */
trait DNASequence {

}


object ACGTSequence {


}

trait DNA2BitEncoding {

  protected def blockIndex(basePos: Long): Int = (basePos >>> 5).toInt

  protected def blockOffset(basePos: Long): Int = (basePos & 0x1FL).toInt

}


/**
 * 2-bit encoded DNA Sequence of A, C, G and T. The maximum size this class can hold is
 * 2G (max of Java array size) * 8 (long byte size) * 8 / 2 (bit) =  64G  (64G characters)
 * @param seq
 * @param numBases
 */
class ACGTSequence(private val seq: Array[Long], val numBases: Long) extends DNASequence with DNA2BitEncoding with CharSequence {
  private val LONG_BYTE_SIZE: Int = 8
  private val MAX_SIZE: Long = (2 * 1024 * 1024 * 1024 * 8 * 8) / 2
  private var hash: Int = 0

  /**
   * length of this sequence
   * @return
   */
  def length = numBases.toInt




  override def toString = {
    val s = new StringBuilder
    for (i <- 0 until numBases)
      s += apply(i).toChar
    s.result
  }

  def apply(index: Long): DNA = {
    // |---   64bit (32 characters) ---|
    // |ACGT|ACGT|ACGT| ...       |ACGT|

    val pos = (index >> 5).toInt // index / 32 (2^5)
    val offset = (index & 0x01FL).toInt // index % 32

    val code = (seq(pos) >>> (offset << 1)) & 0x03
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
        h += lastBlock * 31L;
      }
      hash = h.toInt
    }
    hash
  }

  override def equals(obj: Any): Boolean = {
    obj match {
      case other: ACGTSequence => {
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

  def complement: ACGTSequence = {
    val c = for (b <- seq) yield ~b
    new ACGTSequence(c, numBases)
  }

  /**
   * Create a reverse string of the this sequence.
   *
   * @return Reverse sequence. The returned sequence is NOT a complement of
   *         the original sequence.
   */
  def reverse: ACGTSequence = {
    val buffer = new ACGTSequenceBuffer(numBases)
    for (i <- 0 until numBases) {
      buffer += this(numBases - i - 1)
    }
    buffer.toACGTSequence
  }

  def reverseComplement: ACGTSequence = {
    reverse.complement
  }


  private def fastCount(v: Long, base: DNA): Long = {
    var r = ~0L
    r &= (if ((base.code & 0x02) == 0) ~v else v) >>> 1
    r &= (if ((base.code & 0x01) == 0) ~v else v)
    r &= 0x5555555555555555L
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

      val ePos = blockIndex(end)
      val eOffset = blockOffset(end)

      var count = 0L
      var numAsinMaskedRegion = 0L
      for (pos <- sPos until ePos) {
        var mask: Long = ~0L
        if (pos == sPos) {
          mask <<= sOffset * 2
          numAsinMaskedRegion += sOffset
        }
        if (pos == ePos - 1) {
          mask &= ~(~0L << (eOffset * 2))
          numAsinMaskedRegion += 32 - eOffset
        }
        // Applying bit mask changes all bases in the masked region to As (code=00)
        val v: Long = seq(pos) & mask
        count += fastCount(v, base)
      }
      if (base == DNA.A)
        count - numAsinMaskedRegion
      else
        count
    }

    base match {
      case N => 0L
      case _ => countACGT
    }

  }

  def fastCountACGT(start: Long, end: Long): Array[Long] = {
    val count = new ArrayBuffer[Long](4)

    val sPos = blockIndex(start)
    val sOffset = blockOffset(start)

    val ePos = blockIndex(end)
    val eOffset = blockOffset(end)

    var numAsinMaskedRegion = 0L
    for (pos <- sPos until ePos) {
      var mask: Long = ~0L
      if (pos == sPos) {
        mask <<= sOffset * 2
        numAsinMaskedRegion += sOffset
      }
      if (pos == ePos - 1) {
        mask &= ~(~0L << (eOffset * 2))
        numAsinMaskedRegion += 32 - eOffset
      }
      // Applying bit mask changes all bases in the masked region to As (code=00)
      val v: Long = seq(pos) & mask
      for (base <- DNA.exceptN) {
        count(base.code) += fastCount(v, base)
      }
    }
    count(DNA.A.code) -= numAsinMaskedRegion
    count.toArray
  }

  def charAt(index: Int) = apply(index).toChar

  def subSequence(start: Int, end: Int) = {

  }
}


object ACGTSequenceBuffer {
  private val LONG_BYTE_SIZE: Int = 8
  // 2G (max of Java array size) * 8 (long byte size) * 8 / 2 (bit) =  64G  (64G characters)
  private val MAX_SIZE: Long = 2L * 1024L * 1024L * 1024L * 8L * 8L / 2L

  private def minArraySize(numBases: Long): Int = {
    val bitSize: Long = numBases * 2L
    val blockBitSize: Long = LONG_BYTE_SIZE * 8L
    val arraySize: Long = ((bitSize + blockBitSize - 1L) / blockBitSize) * 2L
    if (arraySize > Integer.MAX_VALUE) {
      throw new IllegalArgumentException("Cannot create ACGTSequece more than %,d characters: %,d".format(MAX_SIZE, numBases))
    }
    arraySize.toInt
  }
}

class ACGTSequenceBuffer(private var seq: ArrayBuffer[Long], var numBases: Long) extends DNASequence with DNA2BitEncoding {

  import ACGTSequenceBuffer._

  private def capacity: Long = seq.length / 2L * 64L

  private var hash = 0 // default to 0

  private def ensureArrayCapacity(newCapacity: Long) {
    if (newCapacity >= capacity) {
      val arraySize = minArraySize(newCapacity)
      val newSeq = new ArrayBuffer[Long](arraySize)
      seq.copyToBuffer(newSeq)
      seq = newSeq
    }
  }

  def +=(base: DNA): Unit = {
    val index = numBases
    numBases += 1
    if (index > capacity) {
      val newCapacity = (index * 3L / 2L) + 64L
      ensureArrayCapacity(newCapacity)
    }
    update(index, base)
  }

  def update(index: Long, base: DNA): Unit = {
    val pos = blockIndex(index)
    val offset = blockOffset(index)
    val shift = offset * 2

    // reset the target base
    seq(pos) &= ~(0x3L << shift)
    seq(pos) |= (base.code & 0x03) << shift
  }

  /**
   * Create an emptry sequence
   */
  def this() = this(new ArrayBuffer[Long](minArraySize(10)), 0L)

  def this(numBases: Long) = this(new ArrayBuffer[Long](numBases), numBases)


  def toACGTSequence: ACGTSequence = new ACGTSequence(seq.toArray, numBases)

  ////   * Create ACGTSequence from the input ACGT(N) sequence
  ////   *
  ////   * @param s
  ////   */
  ////  public ACGTSequence(CharSequence s) {
  ////    this(countNonWhiteSpaces(s));
  ////
  ////    int index = 0;
  ////    for (int i = 0; i < s.length(); ++i) {
  ////      char ch = s.charAt(i);
  ////      if (ch == ' ')
  ////        continue; // skip white space
  ////      byte code = ACGT.to3bitCode(ch);
  ////      set(index++, code);
  ////    }
  ////  }
  //
  //  def replaceN_withA : ACGTSequence =  {
  //    ACGTSequence newSeq = new ACGTSequence(this);
  //    for (int i = 0; i < this.length(); ++i) {
  //      if (this.getACGT(i) == ACGT.N) {
  //        newSeq.set(i, ACGT.A);
  //      }
  //    }
  //    return newSeq;
  //  }
  //
  //  private static long countNonWhiteSpaces(CharSequence s) {
  //    int count = 0;
  //    for (int i = 0; i < s.length(); ++i) {
  //      if (s.charAt(i) != ' ')
  //        count++;
  //    }
  //    return count;
  //  }
  //
  //  /**
  //   * Create a sequence that can hold the given number of bases
  //   *
  //   * @param numBases
  //   */
  //  public ACGTSequence(long numBases) {
  //    this.numBases = numBases;
  //
  //    ensureArrayCapacity(numBases);
  //  }
  //
  //
  //
  //  private ACGTSequence(long[] rawSeq, long numBases) {
  //    this.seq = rawSeq;
  //    this.numBases = numBases;
  //  }
  //
  //  public ACGT getACGT(long index) {
  //    return ACGT.decode((byte) lookup(index));
  //  }
  //
  //  @Override
  //  public long lookup(long index) {
  //    int pos = (int) (index >> 6);
  //    int offset = (int) (index & 0x03FL);
  //    int shift = 62 - ((int) (index & 0x1FL) << 1);
  //
  //    long nFlag = seq[pos * 3] & (1L << (63 - offset));
  //    int code = (int) (seq[pos * 3 + (offset >> 5) + 1] >>> shift) & 0x03;
  //    return nFlag == 0 ? code : 4;
  //  }
  //
  //  public void set(long index, ACGT ch) {
  //    set(index, ch.code);
  //  }
  //
  //
  //
  //  @Override
  //  public int hashCode() {
  //    if (hash != 0)
  //      return hash;
  //    int numFilledBlocks = (int) (numBases / 64L * 3L);
  //    long h = numBases * 31L;
  //    int pos = 0;
  //    for (; pos < numFilledBlocks; ++pos) {
  //      h += seq[pos] * 31L;
  //    }
  //    int offset = (int) (numBases % 64L);
  //    if (offset > 0) {
  //      h += (seq[pos] & (~0L << 64 - offset)) * 31L;
  //      h += (seq[pos + 1] & (offset < 32 ? ~0L << (32 - offset) * 2 : ~0L)) * 31L;
  //      h += (seq[pos + 2] & (offset < 32 ? 0L : ~0L << (64 - offset) * 2)) * 31L;
  //    }
  //    hash = (int) h;
  //    return hash;
  //  }
  //
  //  @Override
  //  public boolean equals(Object obj) {
  //    if (!(obj instanceof ACGTSequence))
  //      return false;
  //
  //    ACGTSequence other = ACGTSequence.class.cast(obj);
  //    if (this.numBases != other.numBases)
  //      return false;
  //
  //    int numFilledBlocks = (int) (numBases / 64L * 3L);
  //    int pos = 0;
  //    for (; pos < numFilledBlocks; ++pos) {
  //      if (this.seq[pos] != other.seq[pos])
  //        return false;
  //    }
  //    int offset = (int) (numBases % 64L);
  //    if (offset > 0) {
  //      long mask[] = new long[3];
  //      mask[0] = ~0L << 64 - offset;
  //      mask[1] = offset < 32 ? ~0L << (32 - offset) * 2 : ~0L;
  //      mask[2] = offset <= 32 ? 0L : ~0L << (64 - offset) * 2;
  //      for (int i = 0; i < mask.length; ++i) {
  //        if ((seq[pos + i] & mask[i]) != (other.seq[pos + i] & mask[i]))
  //          return false;
  //      }
  //    }
  //
  //    return true;
  //  }
  //
  //  /**
  //   * Extract and create a clone of the subsequence of the range [start, end)
  //   *
  //   * @param start
  //   * @param end
  //   * @return
  //   */
  //  public ACGTSequence subString(long start, long end) {
  //    if (start > end)
  //      throw new IllegalArgumentException(String.format("invalid range [%d, %d)", start, end));
  //    final long len = end - start;
  //    int minArraySize = minArraySize(len);
  //    ACGTSequence ss = new ACGTSequence(len);
  //    long[] dest = ss.seq;
  //    Arrays.fill(dest, 0L);
  //
  //    for (long i = 0; i < len;) {
  //      int sPos = (int) ((start + i) >> 6);
  //      int sOffset = (int) ((start + i) & 0x3FL);
  //      int dPos = (int) (i >> 6);
  //      int dOffset = (int) (i & 0x3FL);
  //
  //      int copyLen = 0;
  //      long n = seq[sPos * 3];
  //      long h = seq[sPos * 3 + 1];
  //      long l = seq[sPos * 3 + 2];
  //      if (sOffset == dOffset) {
  //        copyLen = 64;
  //      }
  //      else if (sOffset < dOffset) {
  //        // right shift
  //        int shiftLen = dOffset - sOffset;
  //        copyLen = 64 - dOffset;
  //        // Copy Ns
  //        n >>>= shiftLen;
  //        // Copy ACGT blocks
  //        if (shiftLen < 32) {
  //          l = (h << (64 - shiftLen * 2)) | (l >>> shiftLen * 2);
  //          h >>>= shiftLen * 2;
  //        }
  //        else {
  //          l = h >>> (shiftLen - 32) * 2;
  //          h = 0L;
  //        }
  //      }
  //      else {
  //        // left shift
  //        int shiftLen = sOffset - dOffset;
  //        copyLen = 64 - sOffset;
  //        // Copy Ns
  //        n <<= shiftLen;
  //        // Copy ACGT blocks
  //        if (shiftLen < 32) {
  //          h = (h << shiftLen * 2) | (l >>> (64 - shiftLen * 2));
  //          l <<= shiftLen * 2;
  //        }
  //        else {
  //          h = l << (shiftLen - 32) * 2;
  //          l = 0L;
  //        }
  //      }
  //      dest[dPos * 3] |= n;
  //      dest[dPos * 3 + 1] |= h;
  //      dest[dPos * 3 + 2] |= l;
  //
  //      i += copyLen;
  //    }
  //
  //    return ss;
  //  }
  //
  //  @Override
  //  public long textSize() {
  //    return numBases;
  //  }
  //
  //  /**
  //   * Create a reverse string of the this sequence. The sentinel is appended as
  //   * the last character of the resulting sequence.
  //   *
  //   * @return Reverse sequence. The returned sequence is NOT a complement of
  //   *         the original sequence.
  //   */
  //  public ACGTSequence reverse() {
  //    ACGTSequence rev = new ACGTSequence(this.numBases);
  //    for (long i = 0; i < numBases; ++i) {
  //      rev.set(i, this.lookup(numBases - i - 1));
  //    }
  //    return rev;
  //  }
  //
  //  /**
  //   * Create a complementary sequence (not reversed)
  //   *
  //   * @return complementary sequence
  //   */
  //  public ACGTSequence complement() {
  //    ACGTSequence c = new ACGTSequence(this.numBases);
  //
  //    int numBlocks = seq.length / 3;
  //    for (int i = 0; i < numBlocks; ++i) {
  //      c.seq[i * 3] = this.seq[i * 3];
  //      c.seq[i * 3 + 1] = ~(this.seq[i * 3 + 1]);
  //      c.seq[i * 3 + 2] = ~(this.seq[i * 3 + 2]);
  //    }
  //    return c;
  //  }
  //
  //  public ACGTSequence reverseComplement() {
  //    ACGTSequence rc = reverse();
  //    int numBlocks = seq.length / 3;
  //    for (int i = 0; i < numBlocks; ++i) {
  //      rc.seq[i * 3 + 1] = ~(rc.seq[i * 3 + 1]);
  //      rc.seq[i * 3 + 2] = ~(rc.seq[i * 3 + 2]);
  //    }
  //    return rc;
  //  }
  //
  //  public static ACGTSequence loadFrom(File f) throws IOException {
  //    DataInputStream d = new DataInputStream(new BufferedInputStream(new FileInputStream(f), 4 * 1024 * 1024));
  //    try {
  //      return ACGTSequence.loadFrom(d);
  //    }
  //    finally {
  //      d.close();
  //    }
  //  }
  //
  //  public static ACGTSequence loadFrom(DataInputStream in) throws IOException {
  //    // The num bases must be always 2
  //
  //    long numBases = in.readLong();
  //    int longArraySize = minArraySize(numBases);
  //    long[] seq = new long[longArraySize];
  //    SnappyInputStream sin = new SnappyInputStream(in);
  //    int readBytes = sin.read(seq);
  //    return new ACGTSequence(seq, numBases);
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
  //
  //  @Override
  //  public String toString() {
  //    StringBuilder b = new StringBuilder();
  //    for (int i = 0; i < numBases; ++i) {
  //      ACGT base = ACGT.decode((byte) lookup(i));
  //      b.append(base);
  //    }
  //    return b.toString();
  //  }
  //
  //  /**
  //   * Count the number of the specified character in the range
  //   *
  //   * @param base
  //   * @param start
  //   * @param end
  //   * @return
  //   */
  //  public long count(ACGT base, long start, long end) {
  //    long count = 0;
  //    for (long i = start; i < end; ++i) {
  //      if ((char) lookup(i) == base.code)
  //        count++;
  //    }
  //    return count;
  //  }
  //
  //  /**
  //   * Count the number of occurrence of the code within the specified range
  //   *
  //   * @param base
  //   * @param start
  //   * @param end
  //     *            (exclusive)
  //   * @return
  //   */
  //  public long fastCount(ACGT base, long start, long end) {
  //    long count = 0;
  //
  //    if (base == ACGT.N) {
  //      // Count N
  //      int sPos = (int) (start >>> 6);
  //      int sOffset = (int) (start & 0x3FL);
  //      int ePos = (int) ((end + 64L - 1L) >>> 6);
  //      for (; sPos < ePos; ++sPos) {
  //        long mask = ~0L;
  //        if (sOffset != 0) {
  //          mask >>>= sOffset;
  //          sOffset = 0;
  //        }
  //        if (sPos == ePos - 1) {
  //          int eOffset = (int) (end & 0x3FL);
  //          long rMask = (eOffset == 0) ? ~0L : ~((1L << (64 - eOffset)) - 1);
  //          mask &= rMask;
  //        }
  //        count += Long.bitCount(seq[sPos * 3] & mask);
  //      }
  //    }
  //    else {
  //      // Count A, C, G, T
  //      int sPos = (int) (start >>> 5);
  //      int sOffset = (int) (start & 0x1FL);
  //      int ePos = (int) ((end + 32L - 1L) >>> 5);
  //
  //      for (; sPos < ePos; ++sPos) {
  //
  //        long mask = ~0L;
  //        if (sOffset != 0) {
  //          mask >>>= sOffset * 2;
  //          sOffset = 0;
  //        }
  //        int bIndex = sPos / 2 * 3;
  //        int block = sPos % 2;
  //        long v = seq[bIndex + 1 + block];
  //        long nFlag = interleave32With0(seq[bIndex] >>> (32 * (1 - block)));
  //        if (sPos == ePos - 1) {
  //          int eOffset = (int) (end & 0x1FL);
  //          long rMask = (eOffset == 0) ? ~0L : ~((1L << (32 - eOffset) * 2) - 1);
  //          mask &= rMask;
  //        }
  //        long r = ~0L;
  //        r &= ((base.code & 0x02) == 0 ? ~v : v) >>> 1;
  //        r &= ((base.code & 0x01) == 0 ? ~v : v);
  //        r &= 0x5555555555555555L;
  //        r &= ~nFlag;
  //        r &= mask;
  //        count += Long.bitCount(r);
  //      }
  //    }
  //
  //    return count;
  //  }
  //
  //  public long[] fastCountACGTN(long start, long end) {
  //
  //    long count[] = new long[5];
  //
  //    // Count A, C, G, T
  //    int sPos = (int) (start >>> 5);
  //    int sOffset = (int) (start & 0x1FL);
  //    int ePos = (int) ((end + 32L - 1L) >>> 5);
  //
  //    for (; sPos < ePos; ++sPos) {
  //
  //      long mask = ~0L;
  //      if (sOffset != 0) {
  //        mask >>>= sOffset * 2;
  //        sOffset = 0;
  //      }
  //      int bIndex = sPos / 2 * 3;
  //      int block = sPos % 2;
  //      long v = seq[bIndex + 1 + block];
  //      long nFlag = interleave32With0(seq[bIndex] >>> (32 * (1 - block)));
  //      if (sPos == ePos - 1) {
  //        int eOffset = (int) (end & 0x1FL);
  //        long rMask = (eOffset == 0) ? ~0L : ~((1L << (32 - eOffset) * 2) - 1);
  //        mask &= rMask;
  //      }
  //
  //      for (ACGT base : ACGT.exceptN) {
  //        long r = ~0L;
  //        r &= ((base.code & 0x02) == 0 ? ~v : v) >>> 1;
  //        r &= ((base.code & 0x01) == 0 ? ~v : v);
  //        r &= 0x5555555555555555L;
  //        r &= ~nFlag;
  //        r &= mask;
  //        count[base.code] += Long.bitCount(r);
  //      }
  //
  //      count[ACGT.N.code] += Long.bitCount(nFlag & mask);
  //    }
  //
  //    return count;
  //  }
  //
  //  static int interleaveWith0(int v) {
  //    v = ((v & 0xFF00) << 8) | (v & 0x00FF);
  //    v = ((v << 4) | v) & 0x0F0F0F0F;
  //    v = ((v << 2) | v) & 0x33333333;
  //    v = ((v << 1) | v) & 0x55555555;
  //    return v;
  //  }
  //
  //  /**
  //   * Interleave low 32bits (in a long value) with 0s. For example, 11110011 (8
  //   * bit value) becomes 0101010100000101 (16 bit value)
  //   *
  //   * @param v
  //   * @return
  //   */
  //  static long interleave32With0(long v) {
  //    v = ((v & 0xFFFF0000L) << 16) | (v & 0x0000FFFFL);// 0000000000000000
  //    v = ((v << 8) | v) & 0x00FF00FF00FF00FFL; // 0000000011111111
  //    v = ((v << 4) | v) & 0x0F0F0F0F0F0F0F0FL; // 00001111
  //    v = ((v << 2) | v) & 0x3333333333333333L; // 0011
  //    v = ((v << 1) | v) & 0x5555555555555555L; // 0101
  //    return v;
  //  }
  //
  //  /**
  //   * Count the number of 1s in the input. See also the Hacker's Delight:
  //   * http://hackers-delight.org.ua/038.htm
  //   *
  //   * @param x
  //   * @return the number of 1-bit in the input x
  //   */
  //  public static long countOneBit(long x) {
  //    x = (x & 0x5555555555555555L) + ((x >>> 1) & 0x5555555555555555L);
  //    x = (x & 0x3333333333333333L) + ((x >>> 2) & 0x3333333333333333L);
  //    x = (x + (x >>> 4)) & 0x0F0F0F0F0F0F0F0FL;
  //    x = x + (x >>> 8);
  //    x = x + (x >>> 16);
  //    x = x + (x >>> 32);
  //    return x & 0x7FL;
  //  }
  //
  //  @Override
  //  public long increment(long i, long val) {
  //    throw new UnsupportedOperationException("update");
  //  }
  //
  //  @Override
  //  public int length() {
  //    return (int) textSize();
  //  }
  //
  //  @Override
  //  public char charAt(int index) {
  //    return getACGT(index).toChar();
  //  }
  //
  //  @Override
  //  public ACGTSequence subSequence(int start, int end) {
  //    return subString(start, end);
  //  }
  //
  //}

}