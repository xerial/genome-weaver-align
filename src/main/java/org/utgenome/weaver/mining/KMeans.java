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
// KMeans.java
// Since: 2010/11/30
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.mining;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Random;

import org.xerial.util.ArrayDeque;
import org.xerial.util.Deque;
import org.xerial.util.IndexedSet;
import org.xerial.util.log.Logger;

/**
 * K-means clustering
 * 
 * @author leo
 * 
 */
public class KMeans<T> {

	private static Logger _logger = Logger.getLogger(KMeans.class);

	/**
	 * Distance metric for K-Means clustering
	 * 
	 * @author leo
	 * 
	 */
	public static interface PointMetric<T> {
		/**
		 * Return the distance between the given two points
		 * 
		 * @param a
		 * @param b
		 * @return |a-b|
		 */
		public double distance(T a, T b);

		/**
		 * Return the center of mass of the inputs
		 * 
		 * @param it
		 *            iterator of the list of input points
		 * @return |a+b+c+... | / N
		 */
		public T centerOfMass(Iterator<T> it);

		/**
		 * Return the squared sum of errors of the points from the given center point. \sum_i| x_i - x |^2
		 * 
		 * @param it
		 * @param center
		 * @return
		 */
		public double squareError(Iterator<T> it, T center);

		/**
		 * Compute the lower bound of the points
		 * 
		 * @param it
		 * @return
		 */
		public T lowerBound(Iterator<T> it);

		/**
		 * Compute the upper bound of the points
		 * 
		 * @param it
		 * @return
		 */
		public T upperBound(Iterator<T> it);

		/**
		 * Move the points to the specified direction to the amount of the given distance
		 * 
		 * @param point
		 * @param directionInRadian
		 * @param dist
		 * @return
		 */
		public T move(T point, double directionInRadian, double dist);

		/**
		 * The size of dimension;
		 * 
		 * @return
		 */
		public int dimSize();
	}

	/**
	 * Holds K-means clustering result
	 * 
	 * @author leo
	 * 
	 */
	public static class ClusterInfo<T> {

		/**
		 * the number of clusters;
		 */
		int K;
		/**
		 * Average of squared distance of each point to its belonging centroids
		 */
		public double averageOfDistance = Double.MAX_VALUE;
		/**
		 * List of centroids
		 */
		public IndexedSet<T> centroid = new IndexedSet<T>();

		/**
		 * Array of cluster IDs of the point p_0, ..., p_{N-1};
		 */
		int[] clusterAssignment;

		public ClusterInfo(int K, int N, List<T> centroids) {
			this.K = K;
			this.clusterAssignment = new int[N];
			// register centroids
			for (T each : centroids)
				this.centroid.add(each);
		}
	}

	public static class Config {
		public int maxIteration = 300;
	}

	private final Config config;
	private final PointMetric<T> metric;

	public KMeans(PointMetric<T> metric) {
		this(new Config(), metric);

	}

	public KMeans(Config config, PointMetric<T> metric) {
		this.config = config;
		this.metric = metric;
	}

	/**
	 * @param K
	 *            number of clusters
	 * @param points
	 *            input data points in the matrix format
	 * @param metric
	 *            clustering metric in which the functions for computing distance and the center of mass are defined
	 * @throws Exception
	 */
	public ClusterInfo<T> execute(int K, List<T> points) throws Exception {
		return execute(K, points, initCentroids(K, points));
	}

	private boolean hasDuplicate(List<T> points, T p) {
		for (T each : points) {
			double dist = metric.distance(each, p);
			if (dist == 0)
				return true;
		}
		return false;
	}

	private boolean hasDuplicate(List<T> points) {
		for (int i = 0; i < points.size(); ++i) {
			T a = points.get(i);
			for (int j = i + 1; j < points.size(); ++j) {
				T b = points.get(j);
				if (metric.distance(a, b) == 0.0)
					return true;
			}
		}
		return false;
	}

	/**
	 * Randomly choose K-centroids from the input data set
	 * 
	 * @param K
	 * @param points
	 * @return
	 */
	protected List<T> initCentroids(int K, List<T> points) {
		// N: the number of data points
		final int N = points.size();

		List<T> centroid = new ArrayList<T>(K);

		Random random = new Random(0);

		for (int i = 0; i < K; ++i) {
			int r = random.nextInt(N);

			// TODO avoid infinite loop
			while (hasDuplicate(centroid, points.get(r))) {
				r = random.nextInt(N);
			}
			centroid.add(points.get(r));
		}

		return centroid;
	}

	/**
	 * @param K
	 *            number of clusters
	 * @param points
	 *            input data points in the matrix format
	 * @param centroids
	 *            initial centroids
	 * @param metric
	 *            clustering metric in which the functions for computing distance and the center of mass are defined
	 * @throws Exception
	 */
	public ClusterInfo<T> execute(int K, List<T> points, List<T> centroids) throws Exception {

		// Initialization Step -------------------
		// N: the number of data points
		final int N = points.size();

		// Iteration Step
		ClusterInfo<T> prevClusters = new ClusterInfo<T>(K, N, centroids);
		int maxIteration = config.maxIteration;
		for (int i = 0; i < maxIteration; ++i) {
			if (_logger.isDebugEnabled())
				_logger.debug(String.format("iteration : %d", i + 1));

			// assign each points to the closest centroid 
			ClusterInfo<T> newClusters = EStep(K, points, centroids);
			if (newClusters.averageOfDistance >= prevClusters.averageOfDistance) {
				// Has no improvement in distance 
				return prevClusters;
			}
			// update the centroids
			centroids = MStep(K, points, newClusters);

			if (hasDuplicate(centroids)) {
				// When duplicate centroids are found, stop iteration; 
				return prevClusters;
			}

			// update the cluster set
			prevClusters = newClusters;
		}

		return prevClusters;
	}

	protected ClusterInfo<T> EStep(int K, List<T> points, List<T> centroids) {

		// N: the number of data points
		final int N = points.size();

		if (K != centroids.size())
			throw new IllegalStateException(String.format("K=%d, but # of centrods is %d", K, centroids.size()));

		ClusterInfo<T> result = new ClusterInfo<T>(K, N, centroids);

		// E-step: find the closest centroids, and assign each point to the closest centroid
		int[] clusterAssignment = result.clusterAssignment;
		{
			for (int i = 0; i < N; ++i) {
				T p = points.get(i);
				double dist = Double.MAX_VALUE;
				int closestCentroidID = -1;
				for (int k = 0; k < K; ++k) {
					double d = metric.distance(p, centroids.get(k));
					if (d < dist) {
						dist = d;
						closestCentroidID = k;
					}
				}
				assert (closestCentroidID != -1);
				clusterAssignment[i] = closestCentroidID;
			}
		}

		// Compute the average squared-distance to the centroids
		{
			double averageOfDistance = 0;
			for (int i = 0; i < N; i++) {
				T p = points.get(i);
				T centroidOfP = result.centroid.getByID(clusterAssignment[i]);
				averageOfDistance += Math.pow(metric.distance(p, centroidOfP), 2);
			}
			averageOfDistance /= N;
			result.averageOfDistance = averageOfDistance;
		}

		return result;
	}

	/**
	 * Returns the list of new centroids
	 * 
	 * @param K
	 * @param points
	 * @param clusterInfo
	 * @param metric
	 * @return
	 */
	protected List<T> MStep(int K, List<T> points, ClusterInfo<T> clusterInfo) {

		// M-step: re-estimate the centroids by taking the center of mass of points
		List<T> newCentroid = new ArrayList<T>(K);
		{
			for (int centroidID : clusterInfo.centroid.getIDSet()) {
				T centerOfMass = metric.centerOfMass(new GroupElementsIterator(points, clusterInfo.clusterAssignment, centroidID));
				newCentroid.add(centerOfMass);
			}
		}

		return newCentroid;
	}

	/**
	 * Iterator for selecting points belonging to a specific centroid
	 * 
	 * @author leo
	 * 
	 */
	private class GroupElementsIterator implements Iterator<T> {

		final List<T> points;
		final int[] belongingClass;
		final int k;

		Deque<T> queue = new ArrayDeque<T>();
		int cursor = 0;

		public GroupElementsIterator(List<T> points, int[] belongingClass, int k) {
			this.points = points;
			this.belongingClass = belongingClass;
			this.k = k;
		}

		public boolean hasNext() {
			if (!queue.isEmpty())
				return true;

			for (; cursor < points.size(); ++cursor) {
				if (belongingClass[cursor] == k) {
					T next = points.get(cursor);
					queue.add(next);
					cursor++;
					return true;
				}
			}
			return false;
		}

		public T next() {
			if (hasNext())
				return queue.pollFirst();
			else
				throw new NoSuchElementException();
		}

		public void remove() {
			throw new UnsupportedOperationException("remove");
		}

	}

}
