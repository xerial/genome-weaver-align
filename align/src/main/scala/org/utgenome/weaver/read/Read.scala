package org.utgenome.weaver.read

import org.utgenome.weaver.align.ACGT
import org.utgenome.weaver.align.ACGTSequence
<<<<<<< HEAD:gw-align/src/main/scala/org/utgenome/weaver/read/Read.scala
import java.io.InputStream
import java.io.File
import java.io.FileInputStream
import java.io.BufferedInputStream
import java.io.Reader
import org.utgenome.format.fastq.FastqReader
=======
import java.io.File
import java.io.Reader
>>>>>>> develop:align/src/main/scala/org/utgenome/weaver/read/Read.scala
import java.io.FileReader
import collection._

/**
 * Base class of DNA sequences (string, 2-bit encoding)
 *
 * @author leo
 *
 */
trait DNASequence {
  def apply(i: Long): ACGT

  def size: Long

  def length: Long = size

  def complement: DNASequence

  def count(ch: ACGT, start: Long, end: Long): Long

  def replaceN_withA: DNASequence
}

/**
<<<<<<< HEAD:gw-align/src/main/scala/org/utgenome/weaver/read/Read.scala
 * 2-bit encoding based DNASequence
=======
 * 2-bit encoding based DNASeq
>>>>>>> develop:align/src/main/scala/org/utgenome/weaver/read/Read.scala
 *
 * @author leo
 *
 */
class CompactDNASequence(val seq: ACGTSequence) extends DNASequence {
  def this(seq: String) = this (new ACGTSequence(seq))

  def apply(i: Long): ACGT = seq.getACGT(i)

  def size = seq.textSize();

  def complement = new CompactDNASequence(seq.complement())

  def count(ch: ACGT, start: Long, end: Long) = seq fastCount(ch, start, end)

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
trait Read {
  def numReadFragments: Int

  def apply(i: Int): Read

}

object Read {
  implicit def stringToDNASequence(x: String): DNASequence = new CompactDNASequence(x)

<<<<<<< HEAD:gw-align/src/main/scala/org/utgenome/weaver/read/Read.scala
  //implicit def convertToDNASequence(x: ACGTSequence): DNASequence = new CompactDNASequence(x)
=======
  //implicit def convertToDNASequence(x: ACGTSeq): DNASeq = new CompactDNASequence(x)
>>>>>>> develop:align/src/main/scala/org/utgenome/weaver/read/Read.scala
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

case class PairedEndRead(first: SingleEnd, second: SingleEnd) extends Read {
  def numReadFragments = 2

  def apply(i: Int): SingleEnd = i match {
    case 0 => first
    case 1 => second
  }
}


trait ReadIterator extends Iterator[Read] {
  protected var current: Option[Read] = None
  protected var finishedReading = false

  protected def consume : Option[Read]

  def hasNext: Boolean = consume.isDefined

  def next: Read = {
    val ret = consume match {
      case Some(e) =>  e
      case None => null
    }
    current = None
    ret
  }

  def blockIterator(blockSize:Int=1024): GroupedIterator[Read] = {
     sliding(blockSize, blockSize)
  }

}

<<<<<<< HEAD:gw-align/src/main/scala/org/utgenome/weaver/read/Read.scala
class FASTQFileReader(in: Reader) extends ReadIterator {

  import Read.stringToDNASequence

  val reader = new FastqReader(in)

  def this(file: File) = {
    this (new FileReader(file))
  }

  override def next : SingleEnd = super[ReadIterator].next.asInstanceOf[SingleEnd]

  protected def consume: Option[Read] = {
    if (!current.isDefined && !finishedReading) {
      current = reader.next match {
        case null => finishedReading = true; reader.close; None
        case e => Some(FASTQRead(e.seqname, e.seq, e.qual))
      }
    }
    current
  }

}

class FASTQPairedEndReader(in1: Reader, in2: Reader) extends ReadIterator {
  val reader1 = new FASTQFileReader(in1)
  val reader2 = new FASTQFileReader(in2)

  def this(file1: File, file2: File) = {
    this (new FileReader(file1), new FileReader(file2))
  }

  protected def consume: Option[Read] = {
    if (!current.isDefined && !finishedReading) {
      val r1 = reader1.next
      val r2 = reader2.next
      current = (r1, r2) match {
        case (a: SingleEnd, b: SingleEnd) => Some(PairedEndRead(a, b))
        case (a: SingleEnd, null) => Some(r1)
        case (null, b: SingleEnd) => Some(r2)
        case _ => finishedReading; None
      }
    }
    current
  }

}




=======
//class FASTQFileReader(in: Reader) extends ReadIterator {
//
//  import Read.stringToDNASequence
//
//  val reader = new FastqReader(in)
//
//  def this(file: File) = {
//    this (new FileReader(file))
//  }
//
//  override def next : SingleEnd = super[ReadIterator].next.asInstanceOf[SingleEnd]
//
//  protected def consume: Option[Read] = {
//    if (!current.isDefined && !finishedReading) {
//      current = reader.next match {
//        case null => finishedReading = true; reader.close; None
//        case e => Some(FASTQRead(e.seqname, e.seq, e.qual))
//      }
//    }
//    current
//  }
//
//}
//
//class FASTQPairedEndReader(in1: Reader, in2: Reader) extends ReadIterator {
//  val reader1 = new FASTQFileReader(in1)
//  val reader2 = new FASTQFileReader(in2)
//
//  def this(file1: File, file2: File) = {
//    this (new FileReader(file1), new FileReader(file2))
//  }
//
//  protected def consume: Option[Read] = {
//    if (!current.isDefined && !finishedReading) {
//      val r1 = reader1.next
//      val r2 = reader2.next
//      current = (r1, r2) match {
//        case (a: SingleEnd, b: SingleEnd) => Some(PairedEndRead(a, b))
//        case (a: SingleEnd, null) => Some(r1)
//        case (null, b: SingleEnd) => Some(r2)
//        case _ => finishedReading; None
//      }
//    }
//    current
//  }
//
//}
//
//
//
//
>>>>>>> develop:align/src/main/scala/org/utgenome/weaver/read/Read.scala
