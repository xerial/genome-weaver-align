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

package org.utgenome.weaver
package read

import org.scalatest.FlatSpec
import org.scalatest.matchers.ShouldMatchers
import Read._
import align.ACGT

//--------------------------------------
//
// ReadTest.scala
// Since: 2012/01/04 11:40
//
//--------------------------------------

/**
 * Created at 2012/01/04 11:40
 * @author leo
 */
class ReadTest extends FlatSpec with ShouldMatchers {


  "read" should "handle acgt" in {

    val s :DNASequence = "ACGTACGT"
    val f = new FASTARead("r1", s)

    f.seq.toString should equal (s.toString)
    for(i <- (1L until s.length)) {
      val l : ACGT = f.seq(i)
      val c : ACGT = s(i)
      l should equal (c)
    }
    

  }

}

