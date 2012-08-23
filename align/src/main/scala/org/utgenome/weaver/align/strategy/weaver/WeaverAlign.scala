/*--------------------------------------------------------------------------
 *  Copyright 2011 Taro L. Saito
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *--------------------------------------------------------------------------*/
//--------------------------------------
// genome-weaver project
//
// WeaverAlign.java
// Since: Dec 3, 2011
//
// $URL$
// $Author$
//--------------------------------------

package org.utgenome.weaver
package align.strategy.weaver

import org.utgenome.weaver.align._
import org.utgenome.weaver.read._
import org.utgenome.weaver.parallel.Reporter
import org.utgenome.weaver.read._
import org.utgenome.weaver.align.strategy._
import org.utgenome.weaver.align.SequenceBoundary.PosOnGenome
import scala.collection.mutable.PriorityQueue
import scala.util.control.Breaks._

class Interval(val start: Int, val end: Int) {
  override def toString = "[%d,%d)".format(start, end)

  def length = end - start
}

trait Pivot {
  val pivot: Int
}

trait HasStrand {
  val strand: Strand
}

class ReadRange(val strand: Strand, start: Int, end: Int) extends Interval(start, end) with HasStrand {
  override def toString = "%s%s".format(strand.symbol, super[Interval].toString())
}

/**
 * Points the position of the current search within the read
 * @author leo
 *
 */
sealed abstract class ReadCursor(val readRange: ReadRange, val cursor: Int) {
  def processedBases: Int

  def currentIndex: Int

  def next: Option[ReadCursor]

  def strand = readRange.strand

  def searchDirection: SearchDirection

  def readLength = readRange.length

  override def toString = "%d%s".format(cursor, readRange)
}

class ForwardCursor(readRange: ReadRange, cursor: Int) extends ReadCursor(readRange, cursor) {
  def processedBases = cursor - readRange.start

  def currentIndex = cursor

  def next = {
    val newCursor = cursor + 1
    if (newCursor >= readRange.end) None else Some(new ForwardCursor(readRange, newCursor))
  }

  def searchDirection = SearchDirection.Forward

  override def toString = "F%s".format(super[ReadCursor].toString)
}

class BackwardCursor(readRange: ReadRange, cursor: Int) extends ReadCursor(readRange, cursor) {
  def processedBases = readRange.end - cursor

  def currentIndex = cursor - 1

  def next = {
    val newCursor = cursor - 1
    if (newCursor <= readRange.start) None else Some(new BackwardCursor(readRange, newCursor))
  }

  def searchDirection = SearchDirection.Backward

  override def toString = "B%s".format(super[ReadCursor].toString)
}

class BidirectionalForwardCursor(readRange: ReadRange with Pivot, cursor: Int) extends ReadCursor(readRange, cursor) {
  def processedBases = cursor - readRange.pivot

  def currentIndex = cursor

  def next = {
    val newCursor = cursor + 1
    if (newCursor < readRange.end)
      Some(new BidirectionalForwardCursor(readRange, newCursor))
    else if (readRange.pivot > readRange.start)
      Some(new BackwardCursor(readRange, readRange.pivot))
    else
      None
  }

  def searchDirection = SearchDirection.BidirectionalForward

  //  def atSwitchBoundary = cursor >= readRange.end - 1
  override def toString = "BF%s/%d".format(super[ReadCursor].toString, readRange.pivot)
}

trait Unmapped

trait Mapped {
  val numMismatches: Int
}

trait UniquelyMapped extends Mapped {
  val chr: String
  val start: Int
  val strand: Strand
  val cigar: CIGAR
}

sealed abstract class Alignment(val read: Read)

case class NoHit(override val read: Read) extends Alignment(read) with Unmapped

case class TooManyNs(override val read: Read) extends Alignment(read) with Unmapped

case class ExactMatch(override val read: Read, chr: String, start: Int, strand: Strand, cigar: CIGAR) extends Alignment(read) with UniquelyMapped {
  val numMismatches = 0
}

case class InexactMatch(override val read: Read, chr: String, start: Int, strand: Strand, cigar: CIGAR, numMismatches: Int) extends Alignment(read) with UniquelyMapped

case class FMIndexHit(override val read: Read, si: SuffixInterval, strand: Strand, numMismatches: Int) extends Alignment(read) with Mapped

/**
 *
 * @param fmIndex
 * @param reference
 * @param config
 */
class WeaverAlign(fmIndex: FMIndexOnGenome, reference: ACGTSequence, config: AlignmentConfig) {

  trait Aligner {
    def align(): Alignment
  }

  def align(read: Read): Alignment = {
    val aligner: Aligner = read match {
      case a@FASTARead(name, seq) => new SingleEndAligner(a)
      case b@FASTQRead(name, seq, qual) => new SingleEndAligner(b)
      case p@PairedEndRead(first, second) => new PairedEndAligner(p)
    }

    aligner.align()
  }

  class ReadBreakPoints(
                         val strand: Strand,
                         val si: SuffixInterval,
                         val breakPoints: BitVector,
                         val numMismatches: Int,
                         val longestMatch: Range,
                         val longestMatchSi: SuffixInterval);

  /**
   * Scan the longest unique match region of the query sequence by using FM-index
   * @param query
   * @param strand
   * @return
   */
  def quickScan(query: DNASequence, strand: Strand): ReadBreakPoints = {
    val qLen: Int = query.size.toInt
    val breakPoint = new BitVector(qLen)

    var numMismatches = 0
    var si = fmIndex.wholeSARange()
    var mark = 0
    var longestMatch: Range = null
    var longestMatchSi: SuffixInterval = null;

    var i = 0
    for (i <- (0 until qLen)) {
      val ch = query(i);
      si = fmIndex.forwardSearch(strand, ch, si);
      if (si.isEmpty) {
        breakPoint.set(i, true);
        numMismatches += 1
        if (longestMatch == null || longestMatch.length() < (i - mark)) {
          longestMatch = new Range(mark, i);
          longestMatchSi = si;
        }
        si = fmIndex.wholeSARange();
        mark = i + 1;
      }
    }
    if (longestMatch == null || longestMatch.length() < (i - mark)) {
      longestMatch = new Range(mark, i);
    }
    return new ReadBreakPoints(strand, si, breakPoint, numMismatches, longestMatch, longestMatchSi)
  }

  /**
   * Converts DNASequence for using Java methods that use ACGTSequence
   */
  implicit def dnaToACGT(x: DNASequence): ACGTSequence = {
    x match {
      case a: CompactDNASequence => a.seq
    }
  }

  /**
   * Defines the search order
   */
  implicit val searchOrder = new Ordering[Search] {
    def compare(o1: Search, o2: Search): Int = {
      var diff = 0
      diff = o1.priority - o2.priority
      if (diff == 0)
        diff = -(o1.score - o2.score);
      if (diff == 0)
        diff = -(o1.readCursor.processedBases - o2.readCursor.processedBases);
      return diff;
    }
  }

  class SearchQueue extends PriorityQueue[Search]()(searchOrder) {
    def +=(a: Option[Search]): this.type = {
      a match {
        case Some(x) => this += x
        case None =>
      }
      this
    }

  }

  class SearchStat(var numFMSearch: Int = 0, var numCutOff: Int = 0, var numFiltered: Int = 0)

  class PairedEndAligner(pe: PairedEndRead) extends Aligner {
    def align(): Alignment = {

      val m1 = new SingleEndAligner(pe.first).align
      val m2 = new SingleEndAligner(pe.second).align

      (m1, m2) match {
        case (m1: Mapped, m2: Mapped) =>
        case (u1: Unmapped, m2: Mapped) =>
        case (u1: Unmapped, u2: Unmapped) =>
        case _ => throw new IllegalStateException
      }

      new NoHit(pe)
    }

  }

  val stat = new SearchStat()
  var queryMask: Array[QueryMask] = Array()
  private val staircaseFilterHolder = StaircaseFilter.newHolder

  class SingleEndAligner(read: SingleEnd) extends Aligner {
    private val m: Int = read.length
    private val k: Int = config.getMaximumEditDistance(m)

    def initialState(quickScan: ReadBreakPoints): Option[Search] = {
      if (quickScan.numMismatches > k)
        None;
      else if (quickScan.longestMatch.start != 0 && quickScan.longestMatch.start < m) {
        // bi-directional search
        Some(new Search(new BidirectionalForwardCursor(new ReadRange(quickScan.strand, 0, m) with Pivot {
          val pivot = quickScan.longestMatch.start
        }, quickScan.longestMatch.start),
          k, quickScan.numMismatches))
      } else {
        // single-direction search
        Some(new Search(new ForwardCursor(new ReadRange(quickScan.strand, 0, m), 0), k, quickScan.numMismatches));
      }
    }

    class SearchGraph {

    }

    // Search state queue ordered by priority criteria defined in searchOrder
    private val queue = new SearchQueue

    def align(): Alignment = {
      // Check whether the read contains too many Ns
      def containsTooManyNs(r: SingleEnd) = {
        r.seq.count(ACGT.N, 0, r.length) > k
      }
      if (containsTooManyNs(read)) return new TooManyNs(read)

      // Forward and reverse strand ACGT sequences
      val q: Array[DNASequence] = Array(read.seq.replaceN_withA, read.seq.complement.replaceN_withA)

      // quick scan for k=0.
      def reportExactMatch(rb: ReadBreakPoints): ExactMatch = {
        val pos: PosOnGenome = fmIndex.toGenomeCoordinate(rb.si.lowerBound, m, rb.strand)
        val cigar = new CIGAR("%dM".format(m))
        new ExactMatch(read, pos.chr, pos.pos, pos.strand, cigar)
      }

      // If an exact match is found, return the results immediately
      val scanF = quickScan(q(0), Strand.FORWARD)
      if (scanF.numMismatches == 0) return reportExactMatch(scanF)
      val scanR = quickScan(q(1), Strand.REVERSE)
      if (scanR.numMismatches == 0) return reportExactMatch(scanR)

      // If we are in the exact match mode, no need to proceed furthermore
      if (k == 0) return new NoHit(read);

      // Add initial states for forward and reverse strands
      queue += initialState(scanF)
      queue += initialState(scanR)

      // Create query masks for bit-wise matching
      queryMask = Array(new QueryMask(q(0)), new QueryMask(q(1)))

      // Select a next base or indel to search
      def selectNextSearch(s: Search): NextSearchStep = {
        val nextBase: ACGT = q(s.strandIndex)(s.readCursor.currentIndex)
        // Prefer the next base
        if (!s.isMarked(nextBase)) {
          return Base(nextBase)
        }
        ACGT.exceptN.find {
          ch => !s.isMarked(ch)
        } match {
          case Some(ch) => Base(ch)
          case None => {
            // Lookup indels

            Done()
          }
        }
        Done()
      }

      // Proceed to a next search step
      def stepNext(s: Search, nextCursor: ReadCursor) = {
        selectNextSearch(s) match {
          case Base(nextBase) => {
            val nextState = s.stepOneBase(nextCursor, nextBase)
          }
          case Deletion(length) =>
          case Insertion(length) =>
          case Split() =>
          case SoftClip() =>
          case Done() =>
        }
      }


      while (!queue.isEmpty) {
        val s: Search = queue.dequeue

        breakable {
          // Transit to next state
          s.readCursor.next match {
            case Some(nextCursor) => stepNext(s, nextCursor)
            case None => {
              // report hit
              reportAlignment(FMIndexHit(read, s.fmCursor.currentSi, s.strand, s.nfa.minK))
              break
            }
          }
        }
      }


      new NoHit(read)
    }


    def reportAlignment(aln: Alignment): Unit = {

    }

    def nextSuffixIntervalSet(strand: Strand, cursor: ReadCursor, si: SiSet, currentBase: ACGT) = {
      def nextSi(siF: SuffixInterval, siB: SuffixInterval) = {
        fmIndex.bidirectionalSearch(strand, siF, siB)
      }

      cursor match {
        case f: ForwardCursor => nextSi(si.getForward(currentBase), si.getBackward(currentBase))
        case b: BackwardCursor => nextSi(si.getForward(currentBase), si.getBackward(currentBase))
        case bf: BidirectionalForwardCursor => nextSi(si.getForward(currentBase), si.getBackward(currentBase))
      }
    }

  }

  sealed abstract class NextSearchStep()

  case class Base(base: ACGT) extends NextSearchStep

  case class Deletion(length: Int) extends NextSearchStep

  case class Insertion(length: Int) extends NextSearchStep

  case class Split() extends NextSearchStep

  case class SoftClip() extends NextSearchStep

  case class Done() extends NextSearchStep

  sealed abstract class SearchState

  class SearchEnd extends SearchState

  class Search(val readCursor: ReadCursor, val fmCursor: FMIndexCursor, val nfa: AlignmentNFA, val priority: Int) extends SearchState {

    val nextNFA: Array[Option[AlignmentNFA]] = new Array[Option[AlignmentNFA]](ACGT.exceptN.length)
    (0 until nextNFA.length).foreach {
      nextNFA(_) = None
    }

    var flag: Int = 0

    def this(readCursor: ReadCursor, k: Int, priority: Int) = {
      this (readCursor, new FMIndexCursor(null, fmIndex.initSet(readCursor.searchDirection)), AlignmentNFA.initialNFA(k), priority)
    }

    private def markPos(mark: NextSearchStep) = {
      mark match {
        case Base(base) => base.code
        case Deletion(len) => 6
        case Insertion(len) => 7
        case Split() => 8
        case SoftClip() => 9
        case Done() => 10
      }
    }

    def isMarked(ch: ACGT): Boolean = {
      isMarked(Base(ch))
    }

    def isMarked(mark: NextSearchStep): Boolean = {
      val pos = markPos(mark)
      (flag & (1 << pos)) != 0
    }

    def mark(mark: NextSearchStep): Unit = {
      val pos = markPos(mark)
      flag |= 1 << pos
    }

    //    def getNFA(ch: ACGT): AlignmentNFA = {
    //      if (!nextNFA(ch.code).isDefined) {
    //        val next = nfa.next(cursor, ch, queryMask(cursor.strand.index), staircaseFilterHolder.getStairCaseFilter(cursor.readLength, nfa.minK))
    //        next match {
    //          case None =>
    //          case Some(x) => nextNFA(ch.code) = x.nfa
    //        }
    //      }
    //      nextNFA(ch.code).get
    //    }

    def score: Int = 0

    def strand = readCursor.strand

    def strandIndex: Int = readCursor.strand.index

    def stepOneBase(nextCursor: ReadCursor, base: ACGT): SearchState = {
      val nextFMCursor = fmCursor.next(nextCursor, base)
      stat.numFMSearch += 1

      // TODO update nfa
      new Search(nextCursor, nextFMCursor, nfa, priority)
    }

  }

  class FMIndexCursor(val currentSi: SuffixInterval, val siSet: SiSet) {

    def next(nextCursor: ReadCursor, base: ACGT): FMIndexCursor = {
      val strand = nextCursor.strand

      val nextSi = siSet.getNext(base)
      val nextSiSet: SiSet = nextCursor match {
        case f: ForwardCursor => new SiSet.ForwardSiSet(forwardSearch(strand, nextSi)._1)
        case b: BackwardCursor => new SiSet.BackwardSiSet(backwardSearch(strand, nextSi))
        case fb: BidirectionalForwardCursor => bidirectionalSearch(strand, nextSi, siSet.getBackward(base))
      }
      new FMIndexCursor(nextSi, nextSiSet)
    }

    def forwardSearch(strand: Strand, siF: SuffixInterval) = {
      val fm: FMIndex = if (strand.isForward()) fmIndex.reverseIndex else fmIndex.forwardIndex
      val occLowerBound: Array[Long] = fm.rankACGTN(siF.lowerBound)
      val occUpperBound: Array[Long] = fm.rankACGTN(siF.upperBound)
      val nextSiF: Array[SuffixInterval] = fmIndex.forwardSearch(fm, occLowerBound, occUpperBound)
      (nextSiF, occLowerBound, occUpperBound)
    }

    def backwardSearch(strand: Strand, siB: SuffixInterval) = {
      val fm: FMIndex = if (strand.isForward) fmIndex.forwardIndex else fmIndex.reverseIndex;
      val occLowerBound: Array[Long] = fm.rankACGTN(siB.lowerBound)
      val occUpperBound: Array[Long] = fm.rankACGTN(siB.upperBound)
      fmIndex.forwardSearch(fm, occLowerBound, occUpperBound)
    }

    def bidirectionalSearch(strand: Strand, siF: SuffixInterval, siB: SuffixInterval) = {
      // forward search
      val f = forwardSearch(strand, siF);
      val nextSiB = fmIndex.backwardSearch(f._1, siB, f._2, f._3)
      new SiSet.BidirectionalSiSet(f._1, nextSiB)
    }

  }

}
