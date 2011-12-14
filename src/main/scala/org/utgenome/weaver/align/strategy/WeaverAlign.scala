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

package org.utgenome.weaver.align.strategy

import org.utgenome.weaver.align._
import org.utgenome.weaver.read._
import org.utgenome.weaver.parallel.Reporter
import org.utgenome.weaver.read._
import scala.collection.mutable.PriorityQueue

object WeaverAlign {

}

class Interval(val start: Int, val end: Int) {

}

sealed trait Alignment
case class NoHit(read: Read) extends Alignment
case class TooManyNs(read: Read) extends Alignment
case class FMIndexHit(read: Read, si: SuffixInterval, strand: Strand, numMismatches: Int) extends Alignment

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
   * Scan the longest unique match region in the query sequence by using FM-index
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

  implicit def dnaToACGT(x: DNASequence): ACGTSequence = {
    x match { case a: CompactDNASequence => a.seq }
  }

  class SingleEndAligner(read: SingleEnd) extends Aligner {
    private val m = read.length
    private val k = config.k

    /**
     * Chain of reads
     */
    class Chain(val read: Read) {

      //    class AlignmentState extends Enumeration {
      //      val Init, Mapped, Unmapped = Value
      //    }

      //class Alignment(range: Interval, state: AlignmentState)

    }

    def align(): Alignment = {
      // Check whether the read contains too many Ns
      def containsTooManyNs(r: SingleEnd) = { r.seq.count(ACGT.N, 0, r.length) > k }
      if (containsTooManyNs(read)) return TooManyNs(read)

      val q: Array[DNASequence] = Array(read.seq.replaceN_withA, read.seq.complement.replaceN_withA)

      // quick scan for k=0 
      {
        val wholeReadRange = new Interval(0, m)
        val scanF = quickScan(q(0), Strand.FORWARD)
        if (scanF.numMismatches == 0) return FMIndexHit(read, scanF.si, scanF.strand, 0)
        val scanR = quickScan(q(1), Strand.REVERSE)
        if (scanR.numMismatches == 0) return FMIndexHit(read, scanR.si, scanR.strand, 0)

        if (k == 0) return NoHit(read);

        def selectInitState(quickScan: ReadBreakPoints): Option[Search] = {
          if (scanF.numMismatches > k)
            return None;
          else if (scanF.longestMatch.start != 0 && scanF.longestMatch.start < m) {
            // bi-directional search
            return Some(new Search(new BidirectionalForwardCursor(quickScan.strand, wholeReadRange, quickScan.longestMatch.start, quickScan.longestMatch.start), quickScan.numMismatches))
          } else {
            // single-direction search
            return Some(new Search(new ForwardCursor(quickScan.strand, wholeReadRange), quickScan.numMismatches));
          }
        }

        selectInitState(scanF)

      }

      val queryMask = Array(new QueryMask(q(0)), new QueryMask(q(1)))

      NoHit(read)
    }
  }
  class Search(val cursor: ReadCursor, val priority: Int) {
    def score: Int = 0
  }

  /**
   * Points the position of the current search within the read
   * @author leo
   *
   */
  abstract class ReadCursor(val strand: Strand, val readRange: Interval, val cursor: Int) {
    def processedBases: Int
    def currentIndex: Int
  }

  class ForwardCursor(strand: Strand, readRange: Interval) extends ReadCursor(strand, readRange, readRange.start) {
    def processedBases = readRange.end - cursor
    def currentIndex = cursor
  }

  class BackwardCrusor(strand: Strand, readRange: Interval, cursor: Int) extends ReadCursor(strand, readRange, readRange.end) {
    def processedBases = cursor - readRange.start
    def currentIndex = cursor - 1
  }

  class BidirectionalForwardCursor(strand: Strand, readRange: Interval, cursor: Int, val pivot: Int) extends ReadCursor(strand, readRange, cursor) {
    def processedBases = cursor - pivot
    def currentIndex = cursor
  }

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

  class SearchQueue extends PriorityQueue[Search];

  class PairedEndAligner(pe: PairedEndRead) extends Aligner {
    def align(): Alignment = {

      NoHit(pe)
    }

  }

}
