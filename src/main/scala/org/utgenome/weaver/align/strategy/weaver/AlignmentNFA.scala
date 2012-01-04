package org.utgenome.weaver.align.strategy.weaver

import scala.collection.mutable.StringBuilder
import org.utgenome.weaver.align.ACGT
import org.utgenome.weaver.align.QueryMask
import org.utgenome.weaver.align.strategy.StaircaseFilter

object AlignmentNFA {

  def initialNFA(numAllowedMismatches: Int): AlignmentNFA = {
    new AlignmentNFA(numAllowedMismatches).activateDiagonalStates()
  }

  def toBinary(v: Long, w: Int): String = {
    val s = new StringBuilder
    (0 until w).foreach { i => s append (if (((v >>> i) & 1L) == 0L) "0" else "1") }
    s.toString()
  }
}

/**
 * Non-deterministic finite automaton holding search states of FM-index.
 *
 * <pre>
 *     q0 q1 q2 ... q_{m-1}
 *  k=0 +--+--+--
 *      |\ |\ |\
 *      | \| \| \  ...
 *  k=1 +--+--+--
 *
 * </pre>
 *
 * The states in the automaton are represented as bit flags. kOffset is used to
 * reduce the size of automaton. When no states become active at some layer k,
 * the layer will not be used in future, so we can remove that layer from the
 * automaton.
 *
 * <h3>transitions</h3>
 * <ul>
 * <li>right:match
 * <li>down:deletion from reference
 * <li>right down: insertion to reference (epsilon transition) or mismatch
 * </ul>
 *
 * @author leo
 *
 */
class AlignmentNFA(
  // bit flags holding states at column [index - (k-kOffset), index + (k-kOffset) + 1], where k is # of allowed mismatches 
  automaton: Array[Long],
  kOffset: Int) {
  def this(numAllowedMismatches: Int) = {
    this(new Array[Long](numAllowedMismatches + 1), 0)
  }
  def activateDiagonalStates() = {
    val k: Int = automaton.length - 1;
    (0 until automaton.length) foreach { i => automaton(i) = 1L << (k + i) }
    this
  }

  def minK = kOffset

  override def toString = toNFAStateString
  def toNFAStateString: String = {
    val w = (this.automaton.length + 1) * 2 + 1;
    val s = new StringBuilder();
    s append "kOffset:%d\n".format(kOffset)

    (0 until automaton.length) foreach { j =>
      s.append(AlignmentNFA.toBinary(automaton(j), w));
      if (j != automaton.length - 1) {
        s.append(" \n");
      }
    }
    s.toString();
  }

  def next(cursor: ReadCursor, ch: ACGT, queryMask: QueryMask, staircaseFilter: StaircaseFilter): Option[NextState] = {
    val height = automaton.length;
    val kr = height - 1;
    val pivot = cursor match {
      case p: Pivot => p.pivot
      case _ => 0
    }
    val qeq = queryMask.getBidirectionalPatternMask64(cursor.searchDirection,
      cursor.currentIndex, pivot, cursor.cursor, ch, kr);
    next(cursor.readRange.length, cursor.processedBases, qeq, staircaseFilter)
  }

  def next(readFragmentLength: Int, progressIndex: Int, qeq: Long, staircaseFilter: StaircaseFilter): Option[NextState] = {
    val height: Int = automaton.length
    val k: Int = kOffset + height - 1
    val kr: Int = k - kOffset

    val prev: Array[Long] = automaton
    val next: Array[Long] = new Array(automaton.length);

    var minKwithMatch = k + 1;
    var minKwithProgress = k + 1;

    //String qeqStr = toBinary(qeq, 10);

    //List<Integer> matchPos = new ArrayList<Integer>();
    // Update the automaton
    // R'_0 = ((R_0 & P[ch]) << 1) & (suffix filter)
    next(0) = (prev(0) & qeq) << 1;
    if (next(0) != 0) {
      minKwithMatch = 0;
      minKwithProgress = 0;
      //matchPos.add(cursor.getNextACGTIndex());
    }

    // TODO apply staircase filter

    next(0) &= staircaseFilter.getStairCaseMask64bit(kOffset, progressIndex - k);
    (1 until height).foreach { i =>
      // R'_{i+1} = ((R_{i+1} & P[ch]) << 1) | R_i | (R_i << 1) | (R'_i << 1)
      next(i) = (prev(i) & qeq) << 1;
      if (minKwithMatch > k && next(i) != 0) {
        minKwithMatch = i;

      }
      next(i) |= prev(i - 1) | (prev(i - 1) << 1) | (next(i - 1) << 1);
      // Apply a suffix filter (staircase mask)
      next(i) &= staircaseFilter.getStairCaseMask64bit(kOffset + i, progressIndex - k);
      if (minKwithProgress > k && (next(i) & (1L << height)) != 0L) {
        minKwithProgress = i;
      }
    }

    // Find a match at query position m
    val m = readFragmentLength
    val mPos = k + m - progressIndex
    if (mPos < 64) {
      for (nm <- 0 until height; if ((next(nm) & (1L << mPos)) != 0L))
        return Some(NextState(new AlignmentNFA(removeLayersFromAutomaton(next, nm), kOffset + nm), true))
    }

    // 
    val minK = Math.min(minKwithMatch, minKwithProgress);
    if (minK < k)
      return Some(NextState(new AlignmentNFA(removeLayersFromAutomaton(next, minK), kOffset + minK), false))

    return None;
  }

  case class NextState(nfa: AlignmentNFA, hasMatch: Boolean)

  def removeLayersFromAutomaton(next: Array[Long], numLayersToRemove: Int): Array[Long] = {
    if (numLayersToRemove == 0)
      return automaton;
    val nextHeight = automaton.length - numLayersToRemove;
    // trim rows with fewer mismatches  
    val trimmed = new Array[Long](nextHeight);
    for (h <- 0 until nextHeight)
      trimmed(h) = next(h + numLayersToRemove) >>> 1
    trimmed
  }

  /*
    public ReadAlignmentNFA nextStateAfterSplit(int k) {
        int height = automaton.length - 1;
        long[] nextAutomaton = new long[height];
        for (int i = 0; i < height; ++i) {
            nextAutomaton[i] = 1L << (height + i - 1);
        }
        return new ReadAlignmentNFA(nextAutomaton, kOffset + 1);
    }


    public NextState nextState(Cursor cursor, ACGT ch, QueryMask queryMask, StaircaseFilter staircaseFilter) {
        final int height = automaton.length;
        final int k = kOffset + height - 1;
        final int kr = k - kOffset;
        final long qeq = queryMask.getBidirectionalPatternMask64(cursor.getSearchDirection(),
                cursor.getNextACGTIndex(), cursor.pivot, cursor.cursor, ch, kr);
        return nextState(qeq, cursor.getProcessedBases(), cursor.getFragmentLength(), ch, queryMask, staircaseFilter);
    }

    public NextState nextState(int nextACGTIndex, int progressIndex, int m, ACGT ch, QueryMask queryMask,
            StaircaseFilter staircaseFilter) {
        final int height = automaton.length;
        final int k = kOffset + height - 1;
        final int kr = k - kOffset;
        final long qeq = queryMask.getBidirectionalPatternMask64(SearchDirection.Forward, nextACGTIndex, 0,
                nextACGTIndex, ch, kr);
        return nextState(qeq, progressIndex, m, ch, queryMask, staircaseFilter);
    }

*/
}
