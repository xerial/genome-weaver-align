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
// ImportBEDTest.java
// Since: 2011/07/19
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.db;

import java.io.File;
import java.io.FileReader;

import org.junit.Test;
import org.utgenome.format.bed.BED2Silk;
import org.utgenome.format.bed.BED2SilkReader;
import org.utgenome.util.TestHelper;
import org.utgenome.weaver.GenomeWeaver;
import org.xerial.lens.SilkLens;
import org.xerial.util.log.Logger;

public class ImportBEDTest
{
    private static Logger _logger = Logger.getLogger(ImportBEDTest.class);

    @Test
    public void loadSampleBED() throws Exception {

        File bed = TestHelper.createTempFileFrom(ImportBEDTest.class, "sample.bed");
        GenomeWeaver.execute(String.format("ImportBED %s", bed));

    }

    public static class Track
    {
        public String name;
        public String description;
    }

    public static class BED
    {
        public void addTrack(Track track) {
            _logger.info(SilkLens.toSilk(track));
        }

        public void addGene(BEDAnnotation g) {
            _logger.info(SilkLens.toSilk(g));
        }
    }

    @Test
    public void bed2silk() throws Exception {
        File bed = TestHelper.createTempFileFrom(ImportBEDTest.class, "multi-track.bed");

        BED2Silk bed2silk = new BED2Silk(bed);
        _logger.info("\n" + bed2silk.toSilk());

        SilkLens.loadSilk(BED.class, new BED2SilkReader(new FileReader(bed)));
        _logger.info("");

    }

}
