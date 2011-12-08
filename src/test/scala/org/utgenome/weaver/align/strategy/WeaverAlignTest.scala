/*--------------------------------------------------------------------------
 *  Copyright 2011 utgenome.org
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
// genome-weaver Project
//
// WeaverAlignTest.java
// Since: 2011/12/08
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.strategy

import org.junit._

import Assert._
import org.utgenome.weaver.align._
import org.utgenome.weaver.align.record._
import org.xerial.util.log.Logger

object WeaverAlignTest {
  val _logger: Logger = Logger.getLogger(classOf[BidirectionalSuffixFilterTest])
  val config: AlignmentConfig = new AlignmentConfig()
  val ref: ACGTSequence = new ACGTSequence("AAGCCTAGTTTCCTTG")
  val fmIndex: FMIndexOnGenome = FMIndexOnGenome.buildFromSequence("seq", ref)

  @BeforeClass
  def setUp() {
    //fmIndex = FMIndexOnGenome.buildFromSequence("seq", ref);
    config.k = 2;
  }
}

class WeaverAlignTest {

  @Test
  def hello = {
    println("hello world")
  }

  val p = WeaverAlignTest

  def align(query: String): AlignmentRecord = {
    return align(new ACGTSequence(query));
  }

  def align(q: ACGTSequence): AlignmentRecord = {
    val aln: WeaverAlign = new WeaverAlign(p.fmIndex, p.ref, p.config);
    val result: List[AlignmentRecord] = aln.align(new SingleEndRead("read", q, null))
    if (result.size == 0)
      return null;
    else
      return result(0);
  }

  @Test
  def forwardExact {
    val a: AlignmentRecord = align("GCCTAGT");
    assertEquals("7M", a.cigar.toString());
    assertEquals(3, a.start);
    assertEquals(Strand.FORWARD, a.strand);
    assertEquals(0, a.numMismatches);
  }

  @Test
  def reverseExact = {
    val a: AlignmentRecord = align(new ACGTSequence("GCCTAGT").reverseComplement());
    assertEquals("7M", a.cigar.toString());
    assertEquals(3, a.start);
    assertEquals(Strand.REVERSE, a.strand);
    assertEquals(0, a.numMismatches);
  }

}
