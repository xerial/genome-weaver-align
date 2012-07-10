package xerial.silk.glens


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
// DNASeq.scalaince: 2012/03/16 11:19
//
//--------------------------------------

/**
 * A base trait for classes representing DNA sequences
 *
 * @author leo
 */
trait DNASeq[Repr <: DNASeq[Repr]]
  extends CharSequence
{
  def apply(index:Long) : DNA
  def slice(start:Long, end:Long) : Repr

  def numBases : Long


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


  def count(base:DNA, start:Long, end:Long) : Long
  /**
   * Count the number of occurrences of A, C, G and T letters within [start, end) at the same time.
   * This method is faster than repeating fastCount for each base.
   * @param start
   * @param end
   * @return
   */
  def count(start:Long, end:Long) : Array[Long]

  def charAt(index:Int) : Char = {
    longToIntCheck
    apply(index).toChar
  }

  /**
   * length of this sequence
   * @return
   */
  def length : Int = {
    longToIntCheck
    numBases.toInt
  }

  def subSequence(start:Int, end:Int) = slice(start, end)

  private def longToIntCheck {
    if (numBases >= Integer.MAX_VALUE)
      sys.error("this method cannot be used when the sequence is larger than 2GB")
  }

}

trait DNASeqLike


trait DNASeqOps[Repr <: DNASeq[Repr]]  { this : DNASeq[Repr] =>

  def foreach[U](f:DNA => U) : Unit = {
    for(i <- 0L until numBases)
      f(apply(i))
  }


  /**
   * Take the complement (not reversed) of this sequence
   * @param bf
   * @tparam B
   * @return
   */
  def complement[B](implicit bf:CanBuildSeq[Repr, B]) : B = {
    val b = bf()
    foreach(base => b += base)
    b.result
  }

  /**
   * Create a reverse string of the this sequence. For example ACGT becomes TGCA
   *
   * @return Reverse sequence. The returned sequence is NOT a complement of
   *         the original sequence.
   */
  def reverse[B](implicit bf:CanBuildSeq[Repr, B]) : B = {
    val b = bf()
    var i = 0L
    while (i < numBases) {
      b += apply(numBases - i - 1)
      i += 1
    }
    b.result
  }

  /**
   * Reverse complement of this sequence.
   * @return
   */
  def reverseComplement[B](implicit bf:CanBuildSeq[Repr, B]) : B = {
    val b = bf()
    var i = 0L
    while (i < numBases) {
      b += apply(numBases - i - 1).complement
      i += 1
    }
    b.result
  }

}

trait DNASeqBuilder[Repr] {

  def +=(base:DNA) : Unit
  def +=(seq:String) : Unit = {
    for (ch <- seq)
      this.+=(DNA(ch))
  }

  def result : Repr

}

trait CanBuildSeq[From, To] {
  /**
   *  Creates a new builder on request of a collection.
   *  @param from  the collection requesting the builder to be created.
   *  @return a builder for collections of type `To` with element type `Elem`.
   *          The collections framework usually arranges things so
   *          that the created builder will build the same kind of collection
   *          as `from`.
   */
  def apply(from: From): DNASeqBuilder[To]

  /** Creates a new builder from scratch.
   *
   *  @return a builder for collections of type `To` with element type `Elem`.
   *  @see scala.collection.breakOut
   */
  def apply(): DNASeqBuilder[To]

}
