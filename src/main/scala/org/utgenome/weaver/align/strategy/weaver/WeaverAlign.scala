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

package org.utgenome.weaver.align.strategy.weaver

import org.utgenome.weaver.align._
import org.utgenome.weaver.read._
import org.utgenome.weaver.parallel.Reporter
import org.utgenome.weaver.read._
import org.utgenome.weaver.align.strategy._
import scala.collection.mutable.PriorityQueue

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

sealed trait Alignment
case class NoHit(read: Read) extends Alignment
case class TooManyNs(read: Read) extends Alignment
case class FMIndexHit(read: Read, si: SuffixInterval, strand: Strand, numMismatches: Int) extends Alignment

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
      case a @ FASTARead(name, seq) => new SingleEndAligner(a)
      case b @ FASTQRead(name, seq, qual) => new SingleEndAligner(b)
      case p @ PairedEndRead(first, second) => new PairedEndAligner(p)
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

    var i = 0;
    for (i <- (0 until qLen)) {
      val ch = query(i);
      si = fmIndex.forwardSearch(strand, ch, si);
      if (si.isEmpty()) {
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
    x match { case a: CompactDNASequence => a.seq }
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
        diff = -(o1.cursor.processedBases - o2.cursor.processedBases);
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

  val stat = new SearchStat
  var queryMask: Array[QueryMask]
  private val staircaseFilterHolder = StaircaseFilter.newHolder

  class SingleEndAligner(read: SingleEnd) extends Aligner {
    private val m: Int = read.length
    private val k: Int = config.getMaximumEditDistances(m)

    def initialState(quickScan: ReadBreakPoints): Option[Search] = {
      if (quickScan.numMismatches > k)
        None;
      else if (quickScan.longestMatch.start != 0 && quickScan.longestMatch.start < m) {
        // bi-directional search
        Some(new Search(new BidirectionalForwardCursor(new ReadRange(quickScan.strand, 0, m) with Pivot { val pivot = quickScan.longestMatch.start }, quickScan.longestMatch.start),
          k, quickScan.numMismatches))
      } else {
        // single-direction search
        Some(new Search(new ForwardCursor(new ReadRange(quickScan.strand, 0, m), 0), k, quickScan.numMismatches));
      }
    }

    // Search state queue ordered by priority criteria defined in searchOrder
    private val queue = new SearchQueue

    def align(): Alignment = {
      // Check whether the read contains too many Ns
      def containsTooManyNs(r: SingleEnd) = { r.seq.count(ACGT.N, 0, r.length) > k }
      if (containsTooManyNs(read)) return TooManyNs(read)

      // Forward and reverse strand ACGT sequences
      val q: Array[DNASequence] = Array(read.seq.replaceN_withA, read.seq.complement.replaceN_withA)

      // quick scan for k=0 
      val scanF = quickScan(q(0), Strand.FORWARD)
      if (scanF.numMismatches == 0) return FMIndexHit(read, scanF.si, scanF.strand, 0)
      val scanR = quickScan(q(1), Strand.REVERSE)
      if (scanR.numMismatches == 0) return FMIndexHit(read, scanR.si, scanR.strand, 0)

      // When some exact match is found
      if (k == 0) return NoHit(read);

      // Add initial states for both strand
      queue += initialState(scanF)
      queue += initialState(scanR)

      queryMask = Array(new QueryMask(q(0)), new QueryMask(q(1)))

      while (!queue.isEmpty) {
        val s: Search = queue.dequeue

        // Select a next base or indel to search
        def selectNextSearch(): NextSearchStep = {
          val nextBase: ACGT = q(s.strandIndex)(s.cursor.currentIndex)
          if (!s.isVisited(nextBase)) {
            return Base(nextBase)
          }
          ACGT.exceptN.find { ch => !s.isVisited(ch) } match {
            case Some(ch) => Base(ch)
            case None => Done()
          }
        }

        // Transit to next state
        def stepNext(nextCursor: ReadCursor) = {
          selectNextSearch() match {
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

        s.cursor.next match {
          case Some(nextCursor) => stepNext(nextCursor)
          case None => { // report hit

          }
        }
      }

      NoHit(read)
    }

    sealed class NextSearchStep
    case class Base(base: ACGT) extends NextSearchStep
    case class Deletion(length: Int) extends NextSearchStep
    case class Insertion(length: Int) extends NextSearchStep
    case class Split extends NextSearchStep
    case class SoftClip extends NextSearchStep
    case class Done extends NextSearchStep

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

  sealed abstract class SearchState
  class SearchEnd extends SearchState

  class Search(val cursor: ReadCursor, val fmCursor: FMIndexCursor, val nfa: AlignmentNFA, val priority: Int) extends SearchState {

    val nextNFA: Array[Option[AlignmentNFA]] = new Array[Option[AlignmentNFA]](ACGT.exceptN.length)
    (0 until nextNFA.length).foreach { nextNFA(_) = None }

    def this(cursor: ReadCursor, k: Int, priority: Int) = {
      this(cursor, new FMIndexCursor(null, fmIndex.initSet(cursor.searchDirection)), AlignmentNFA.initialNFA(k), priority)
    }

    def isVisited(ch: ACGT): Boolean = {
      nextNFA(ch.code).isDefined
    }

    def getNFA(ch: ACGT): AlignmentNFA = {
      if (!nextNFA(ch.code).isDefined) {
        val next = nfa.next(cursor, ch, queryMask(cursor.strand.index), staircaseFilterHolder.getStairCaseFilter(cursor.readLength, nfa.minK))
        next match {
          case None =>
          case Some(x) => nextNFA(ch.code) = x.nfa
        }
      }
      nextNFA(ch.code).get
    }

    def score: Int = 0
    def strand = cursor.strand
    def strandIndex: Int = cursor.strand.index

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

  class PairedEndAligner(pe: PairedEndRead) extends Aligner {
    def align(): Alignment = {

      NoHit(pe)
    }

  }

}
