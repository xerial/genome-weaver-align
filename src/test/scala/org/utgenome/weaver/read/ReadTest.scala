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
import Read._
import align.ACGT
import org.xerial.util.FileResource
import org.scalatest.matchers.{MustMatchers, ShouldMatchers}

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
class ReadTest extends FlatSpec with ShouldMatchers with MustMatchers with Logger {

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

  "fastq reader" should "read FASTQ files" in {

    val r = new FASTQFileReader(FileResource.open(classOf[ReadTest], "sample.fastq"))
    val read = r.toArray

    read.length must be (3)
    for(i <- 0 until  read.length) {
      read(i).numReadFragments must equal (1)
    }
  }

  "fastq reader" should "read data in block-wise manner" in {

    val r = new FASTQFileReader(FileResource.open(classOf[ReadTest], "sample.fastq"))
    val block = r.blockIterator(2).toArray

    block.size must be (2)
    block(0).size must be (2)
    block(1).size must be (1)
    val numReads = block.flatMap(_.iterator).size
    numReads must be (3)

    debug(block(0).toString)
    debug(block(1).toString)
  }


}

