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
import org.utgenome.weaver.align.record.Read
import org.utgenome.weaver.parallel.Reporter
import org.utgenome.weaver.align.record.AlignmentRecord

object WeaverAlign {

}

class Interval(start: Int, end: Int) {

}

class WeaverAlign(fmIndex: FMIndexOnGenome, reference: ACGTSequence, config: AlignmentConfig) {

  /**
   * Chain of reads
   */
  class Chain(val read: Read) {

    class AlignmentState extends Enumeration {
      val Init, Mapped, Unmapped = Value
    }

    class Alignment(range: Interval, state: AlignmentState)

  }

  def align(input: Read): List[AlignmentRecord] = {
    val read = for (i <- (0 until input.getNumReadFragment())) yield { input.getRead(i) }

    return List()
  }

}
