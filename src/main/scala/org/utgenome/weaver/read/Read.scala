package org.utgenome.weaver.read
import org.utgenome.weaver.align.ACGT
import org.utgenome.weaver.align.ACGTSequence
import java.io.InputStream
import java.io.File
import java.io.FileInputStream
import java.io.BufferedInputStream
import java.io.Reader
import org.utgenome.format.fastq.FastqReader
import java.io.FileReader

/**
 * Base class of DNA sequences (string, 2-bit encoding)
 *
 * @author leo
 *
 */
trait DNASequence {
  def apply(i: Long): ACGT
  def size: Long
  def complement: DNASequence
  def count(ch: ACGT, start: Long, end: Long): Long

  def replaceN_withA: DNASequence
}

/**
 * 2-bit encoding based DNASequence
 *
 * @author leo
 *
 */
class CompactDNASequence(val seq: ACGTSequence) extends DNASequence {
  def this(seq: String) = this(new ACGTSequence(seq))
  def apply(i: Long): ACGT = seq.getACGT(i)
  def size = seq.textSize();

  def complement = new CompactDNASequence(seq.complement())
  def count(ch: ACGT, start: Long, end: Long) = seq fastCount (ch, start, end)

  def replaceN_withA = {
    new CompactDNASequence(seq.replaceN_withA());
  }

  override def toString = seq.toString
}

/**
 * Read Sequence
 *
 * @author leo
 *
 */
trait Read extends ReadConverter {
  def numReadFragments: Int
  def apply(i: Int): Read

}

trait ReadConverter {
  implicit def stringToDNASequence(x: String): DNASequence = new CompactDNASequence(x)
  //implicit def convert(x: ACGTSequence): DNASequence = new CompactDNASequence(x)
}

trait SingleEnd extends Read {
  val seq: DNASequence
  val length: Int = seq.size.toInt
}

case class FASTARead(name: String, seq: DNASequence) extends SingleEnd {
  def numReadFragments = 1
  def apply(i: Int) = this
}
case class FASTQRead(name: String, seq: DNASequence, val qual: String) extends SingleEnd {
  def numReadFragments = 1
  def apply(i: Int) = this
}
case class PairedEndRead(val first: SingleEnd, val second: SingleEnd) extends Read {
  def numReadFragments = 2
  def apply(i: Int): SingleEnd = i match {
    case 0 => first
    case 1 => second
  }
}

trait ReadStream extends ReadConverter {
  def next: Option[Read]
}

class FASTQStream(in: Reader) extends ReadStream {
  val reader = new FastqReader(in)
  def this(file: File) = {
    this(new FileReader(file))
  }
  def next: Option[SingleEnd] = {
    reader.next match {
      case null => None
      case e => Some(FASTQRead(e.seqname, e.seq, e.qual))
    }
  }
}

class FASTQPairedEndStream(in1: Reader, in2: Reader) extends ReadStream {
  val reader1 = new FASTQStream(in1)
  val reader2 = new FASTQStream(in2)

  def this(file1: File, file2: File) = {
    this(new FileReader(file1), new FileReader(file2))
  }

  def next = {
    val r1 = reader1.next
    val r2 = reader2.next
    (r1, r2) match {
      case (Some(a), Some(b)) => Some(PairedEndRead(a, b))
      case (Some(x), None) => r1
      case (None, Some(y)) => r2
      case _ => None
    }
  }
}

