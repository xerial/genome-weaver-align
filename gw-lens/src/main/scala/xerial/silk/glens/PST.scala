package xerial.silk.glens

import annotation.tailrec
import collection.mutable.Stack
import collection.immutable.{TreeSet, SortedSet}

//--------------------------------------
//
// PST.scala
// Since: 2012/07/10 12:42
//
//--------------------------------------

class Entry[A](val elem: A, val x: Int, val y: Int)

abstract class Node[+A] {
  def isEmpty : Boolean
  def foreach[U](f:(Node[A]) => U) : U
}

/**
 * Persistent priority search tree implementation.
 *
 *
 * @author leo
 */
class PST[+A](protected val root: Node[A], val size: Int) {

  case object Empty extends Node[Nothing] {
    def isEmpty = true
    def foreach[U](f: (Node[Nothing]) => U) {}
  }

  case class Tree(entry: Entry[A], ySplit: Int, left: Node[A], right: Node[A]) extends Node[A] {
    def isEmpty = false
    def foreach[U](f: (Node[A]) => U) {
      f(this)
    }
  }

  def insert(v: A, start: Int, end: Int): PST[A] = {

    def insert(target: Node[A], n: Entry[A], y_lower: Int, y_upper: Int): Node[A] = {
      target match {
        case Empty => Tree(n, (y_lower + y_upper) / 2, Empty, Empty)
        case Tree(c: Entry[A], ySplit, left, right) => {
          // c: current node
          def insertNode(p: Entry[A], e: Entry[A]): Node[A] = {
            if (e.y < ySplit)
              Tree(p, ySplit, insert(left, e, y_lower, ySplit), right)
            else
              Tree(p, ySplit, left, insert(right, e, ySplit, y_upper))
          }

          // swap the entries so that low x value node becomes the parent
          if (n.x < c.x)
            insertNode(n, c) // push down current node
          else
            insertNode(c, n) // insert under the current node
        }
      }
    }

    val newRoot = insert(root, new Entry(v, start, end), 0, Int.MaxValue)
    new PST(newRoot, size + 1)
  }

  def remove(v: A, x: Int, y: Int): PST[A] = {

    var removed = false

    def find(n: Node[A], y1: Int, y2: Int): Node[A] = {
      n match {
        case Empty => Empty
        case t @ Tree(c: Entry[A], ySplit, left, right) =>
          if (v == c.elem) {
            removed = true
            remove(t)
          }
          else {
            if (y < ySplit)
              Tree(c, ySplit, find(left, y1, ySplit), right)
            else
              Tree(c, ySplit, left, find(right, ySplit, y2))
          }
      }
    }

    def remove(t: Tree): Node[A] = {
       (t.left, t.right) match {
        case (left @ Tree(l, _, _, _), right @ Tree(r, _, _, _)) =>
          if(l.x < r.x)
            Tree(l, t.ySplit, remove(left), right) // pull up left
          else
            Tree(r, t.ySplit, left, remove(right)) // pull up right
        case (left @ Tree(l, _, _, _), Empty) => Tree(l, t.ySplit, remove(left), Empty) // pull up left
        case (Empty, right @ Tree(r, _, _, _)) => Tree(r, t.ySplit, Empty, remove(right)) // pull up right
        case (Empty, Empty) => Empty
      }
    }

    val newRoot = find(root, 0, Int.MaxValue)
    if(removed)
      new PST(newRoot, size-1)
    else
      this
  }


  def overlapQuery(start: Int, end: Int): Seq[A] = rangeQuery(end - 1, start + 1, Int.MaxValue)

  def rangeQuery(upperX: Int, ymin: Int, ymax: Int): Seq[A] = {
    val b = Seq.newBuilder[A]
    private class Context(val n: Node[A], val y1: Int, val y2: Int)

    val stack = new Stack[Context]
    if(!root.isEmpty)
      stack.push(new Context(root, ymin, ymax))

    while (!stack.isEmpty) {
      val c: Context = stack.pop
      for(t @ Tree(e, ySplit, left, right) <- c.n) {
        if (e.x <= upperX) {
          if (c.y1 <= e.y && e.y <= c.y2)
            b += e.elem

          if(!left.isEmpty)
            stack.push(new Context(left, c.y1, ySplit))
          if(!right.isEmpty)
            stack.push(new Context(right, ySplit, c.y2))
        }
      }
    }

    b.result
  }


}