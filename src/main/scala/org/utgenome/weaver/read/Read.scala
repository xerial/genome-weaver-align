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
}

/**
 * String-based DNASequence
 *
 * @author leo
 *
 */
class StringDNASequence(seq: String) extends DNASequence {
  def apply(i: Long): ACGT = {
    ACGT.encode(seq.charAt(i.toInt))
  }
  def size = seq.size
}

/**
 * 2-bit encoding based DNASequence
 *
 * @author leo
 *
 */
class CompactSequence(seq: ACGTSequence) extends DNASequence {
  def apply(i: Long): ACGT = seq.getACGT(i)
  def size = seq.textSize();
}

/**
 * Read Sequence
 *
 * @author leo
 *
 */
trait Read {

}

case class FASTARead(val name: String, val seq: DNASequence) extends Read
case class FASTQRead(val name: String, val seq: DNASequence, val qual: String) extends Read
case class PairedEndRead(val first: Read, val second: Read) extends Read

trait ReadStream {
  def next: Option[Read]
}

class FASTQStream(in: Reader) extends ReadStream {
  val reader = new FastqReader(in)
  def this(file: File) = {
    this(new FileReader(file))
  }
  def next = {
    val e = reader.next()
    if (e == null)
      None
    else
      Some(FASTQRead(e.seqname, new StringDNASequence(e.seq), e.qual))
  }
}

class FASTQPairedEndStream(in1: Reader, in2: Reader) extends ReadStream {
  val reader1 = new FASTQStream(in1)
  val reader2 = new FASTQStream(in2)

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

