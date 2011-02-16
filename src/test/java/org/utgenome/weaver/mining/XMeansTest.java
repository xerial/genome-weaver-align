/*--------------------------------------------------------------------------
 *  Copyright 2008 utgenome.org
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
// tss-toolkit Project
//
// XMeansTest.java
// Since: 2010/11/30
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.mining;

import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.List;

import org.junit.Test;
import org.utgenome.weaver.mining.XMeans.ClusterInfo;
import org.xerial.util.FileResource;
import org.xerial.util.log.Logger;

public class XMeansTest
{

    private static Logger _logger = Logger.getLogger(XMeansTest.class);

    @Test
    public void faithful() throws Exception {
        List<Point2D> input = new ArrayList<Point2D>();

        // read the faithful data set
        BufferedReader faithful = FileResource.open(XMeansTest.class, "faithful.txt");
        try {
            for (String line; (line = faithful.readLine()) != null;) {
                String[] c = line.split("\\s+");
                input.add(new Point2D.Double(Double.parseDouble(c[0]), Double.parseDouble(c[1])));
            }

            XMeans<Point2D> xmeans = new XMeans<Point2D>(new EuclidDistanceMetric2D());
            ClusterInfo<Point2D> result = xmeans.execute(input, 1000);

            _logger.info("number of clusters: " + result.K);
            for (int i = 0; i < input.size(); ++i) {
                Point2D p = input.get(i);
                int clusterID = result.clusterAssignment[i];
                System.out.println(String.format("%f\t%f\t%d", p.getX(), p.getY(), clusterID));
            }

            for (Point2D p : result.centroid) {
                System.out.println(String.format("%f\t%f\t%d", p.getX(), p.getY(), -1));
            }
        }
        finally {
            if (faithful != null)
                faithful.close();
        }

    }

}
