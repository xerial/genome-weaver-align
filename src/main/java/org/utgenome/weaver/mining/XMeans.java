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
// XMeans.java
// Since: 2010/11/30
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.mining;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.utgenome.weaver.mining.KMeans.PointMetric;
import org.xerial.util.IndexedSet;

/**
 * X-Means clustering
 * 
 * TODO: normalization of each point
 * 
 * @author leo
 * 
 */
public class XMeans<T>
{

    private final PointMetric<T> metric;

    public XMeans(PointMetric<T> metric) {
        this.metric = metric;
    }

    public static class ClusterInfo<T>
    {
        /**
         * The number of clusters
         */
        public int           K;
        /**
         * List of centroids
         */
        public IndexedSet<T> centroid = new IndexedSet<T>();
        /**
         * Array of the cluster IDs assigned for the points p_0, ... , p_{N-1}
         */
        public int[]         clusterAssignment;

        public double        BIC      = 0;

        public ClusterInfo(int N, T centroid) {
            this.K = 1;
            this.clusterAssignment = new int[N];
            this.centroid.add(centroid);
            for (int i = 0; i < N; i++) {
                clusterAssignment[i] = 0;
            }
        }

        public ClusterInfo(int N, List<T> centroid) {
            this.K = centroid.size();
            this.clusterAssignment = new int[N];
            this.centroid.addAll(centroid);
        }

        public ClusterInfo(int N, KMeans.ClusterInfo<T> other) {
            this.K = other.K;
            this.centroid = other.centroid;
            clusterAssignment = new int[N];
            System.arraycopy(other.clusterAssignment, 0, clusterAssignment, 0, other.clusterAssignment.length);
        }

    }

    /**
     * Compute Bayesian Information Criteria (BIC) of the clusters
     * 
     * @return BIC value
     */
    protected double computeBIC(ClusterInfo<T> cluster, List<T> points) {

        final double R = points.size();
        final double K = cluster.K;

        if (R <= cluster.K) {
            return Double.MIN_VALUE;
        }

        // Compute the sigma of the points under the identical spherical Gaussian assumption
        double sumOfSquaredError = 0.0;
        for (int i = 0; i < points.size(); ++i) {
            int centroidID = cluster.clusterAssignment[i];
            T centroid = cluster.centroid.getByID(centroidID);
            sumOfSquaredError += Math.pow(metric.distance(points.get(i), centroid), 2);
        }
        double sigmaSquare = sumOfSquaredError / (R - K);

        double BIC = 0.0;
        final int M = metric.dimSize();
        for (int k = 0; k < K; ++k) {
            // count the number of elements in the cluster k
            double R_n = 0;
            for (int i = 0; i < points.size(); ++i) {
                if (cluster.clusterAssignment[i] == k) {
                    R_n++;
                }
            }

            // Compute the likelihood of the cluster k
            double p1 = -((R_n / 2.0) * Math.log(2.0 * Math.PI));
            double p2 = -((R_n * M) / 2.0) * Math.log(sigmaSquare);
            double p3 = -(R_n - K) / 2.0;
            double p4 = R_n * Math.log(R_n);
            double p5 = -R_n * Math.log(R);
            double likelihoodOfTheCluster = p1 + p2 + p3 + p4 + p5;

            int numberOfFreeParameters = (int) ((K - 1) + M * K + 1);
            BIC += likelihoodOfTheCluster - (numberOfFreeParameters / 2.0) * Math.log(R_n);
        }

        if (Double.isNaN(BIC))
            return Double.MIN_VALUE;
        return BIC;
    }

    /**
     * @param points
     * @param maxK
     *            maximum number of the cluster
     * @return
     * @throws Exception
     */
    public ClusterInfo<T> execute(List<T> points, int maxK) throws Exception {

        // Start with a single centroid
        T centroid = metric.centerOfMass(points.iterator());
        ClusterInfo<T> cluster = new ClusterInfo<T>(points.size(), centroid);
        cluster.BIC = computeBIC(cluster, points);
        return iteration(cluster, points, maxK);
    }

    protected ClusterInfo<T> iteration(ClusterInfo<T> currentCluster, List<T> points, int maxK) throws Exception {

        for (; currentCluster.K < maxK;) {
            List<T> nextCentroids = findCentroids(currentCluster.K, points, currentCluster);
            ClusterInfo<T> nextCluster = EStep(nextCentroids.size(), points, nextCentroids);

            if (currentCluster.BIC >= nextCluster.BIC)
                break;

            currentCluster = nextCluster;
        }
        return currentCluster;
    }

    protected ClusterInfo<T> EStep(int K, List<T> points, List<T> centroids) throws Exception {
        // Clustering Phase
        KMeans<T> kmeans = new KMeans<T>(metric);
        ClusterInfo<T> cluster = new ClusterInfo<T>(points.size(), kmeans.execute(K, points, centroids));
        cluster.BIC = computeBIC(cluster, points);
        return cluster;
    }

    protected List<T> findCentroids(int K, List<T> points, ClusterInfo<T> cluster) throws Exception {

        // Split Phase
        List<T> nextCentroids = new ArrayList<T>();
        for (int k = 0; k < K; ++k) {
            // Extract the input points that belong to the cluster k 
            List<T> pointsInTheCluster = null;
            if (K != 1) {
                pointsInTheCluster = new ArrayList<T>();
                for (int i = 0; i < points.size(); ++i) {
                    if (cluster.clusterAssignment[i] == k)
                        pointsInTheCluster.add(points.get(i));
                }
            }
            else
                pointsInTheCluster = points;

            // Compute the current BIC
            T centroid = cluster.centroid.getByID(k);

            if (pointsInTheCluster.size() == 1) {
                // When it is a single point cluster
                nextCentroids.add(centroid);
                continue;
            }
            else if (pointsInTheCluster.size() == 0) {
                // When no cluster is found
                continue;
            }

            ClusterInfo<T> currentCluster = new ClusterInfo<T>(pointsInTheCluster.size(), centroid);
            double currentBIC;
            if (K != 1)
                currentBIC = computeBIC(currentCluster, pointsInTheCluster);
            else
                currentBIC = cluster.BIC;

            // Try to improve the cluster structure by splitting the cluster into two 
            ClusterInfo<T> newCluster = split(centroid, pointsInTheCluster);
            if (newCluster == null) {
                nextCentroids.add(centroid);
                continue;
            }

            // Compute the BIC of the split cluster
            double newBIC = computeBIC(newCluster, pointsInTheCluster);

            if (newBIC > currentBIC) {
                // The new cluster has better BIC 
                nextCentroids.addAll(newCluster.centroid);
            }
            else {
                // Using the current centroid as is.
                nextCentroids.add(centroid);
            }
        }

        // remove duplicates
        return removeDuplicates(nextCentroids);
    }

    private List<T> removeDuplicates(List<T> list) {

        List<T> noDup = new ArrayList<T>();
        for (int i = 0; i < list.size(); ++i) {
            T a = list.get(i);
            boolean hasDuplicate = false;
            for (int j = 0; j < list.size(); ++j) {
                if (i == j)
                    continue;
                T b = list.get(j);
                if (metric.distance(a, b) == 0.0) {
                    hasDuplicate = true;
                    break;
                }
            }
            if (!hasDuplicate)
                noDup.add(a);
        }

        return noDup;
    }

    protected ClusterInfo<T> split(T centroid, List<T> points) throws Exception {

        // Compute the range of the cluster
        T lowerBound = metric.lowerBound(points.iterator());
        T upperBound = metric.upperBound(points.iterator());

        // Determine the new centroids by choosing points distant from the centroid at the random direction  
        double diameter = metric.distance(lowerBound, upperBound);
        if (Double.isInfinite(diameter))
            diameter = Double.MAX_VALUE / 2.0;

        if (diameter == 0.0) {
            return null;
        }

        Random rand = new Random(0); // uses a fixed seed for stabilizing the clustering results
        double direction = rand.nextDouble() * Math.PI;
        List<T> newCentroid = new ArrayList<T>();
        newCentroid.add(metric.move(centroid, direction, diameter / 2.0));
        newCentroid.add(metric.move(centroid, direction, -diameter / 2.0));

        // Compute the K-means for the splits
        KMeans<T> kmeans = new KMeans<T>(metric);
        ClusterInfo<T> cluster = new ClusterInfo<T>(points.size(), kmeans.execute(newCentroid.size(), points,
                newCentroid));

        return cluster;
    }

}
