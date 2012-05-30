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
// SOLiDColorTest.java
// Since: 2011/10/12
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import static org.junit.Assert.*;

import java.util.ArrayList;

import org.junit.Test;
import org.xerial.util.CollectionUtil;
import org.xerial.util.Functor;
import org.xerial.util.StringUtil;
import org.xerial.util.log.Logger;

public class SOLiDColorTest
{
    private static Logger _logger = Logger.getLogger(SOLiDColorTest.class);

    @Test
    public void genColorToCuodeTable() throws Exception {

        _logger.debug(
                "{%s}",
                StringUtil.join(
                        CollectionUtil.collect(SOLiDColor.genCharToColorCodeTable(), new Functor<SOLiDColor, String>() {
                            @Override
                            public String apply(SOLiDColor input) {
                                return input.name();
                            }
                        }), ", "));

    }
}
