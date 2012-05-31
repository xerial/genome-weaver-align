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
// MergeSort.java
// Since: 2011/04/08
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.dc3;

import java.util.Comparator;

public class MergeSort
{

    static <T> void swap(T[] array, int x, int y) {
        T tmp = array[x];
        array[x] = array[y];
        array[y] = tmp;
    }

    public static <T> void mergeSort(T[] input, Comparator<T> comparator) {
        T[] dest = input.clone();
        mergeSort(dest, input, 0, input.length, 0, comparator);
    }

    private static int INSERTIONSORT_THRESHOLD = 7;

    static <T> void mergeSort(T[] dest, T[] src, int low, int high, int off, Comparator<T> comparator) {

        int length = high - low;
        // Insertion sort on smallest arrays
        if (length < INSERTIONSORT_THRESHOLD) {
            for (int i = low; i < high; i++)
                for (int j = i; j > low && comparator.compare(dest[j - 1], dest[j]) > 0; j--) {
                    swap(dest, j - 1, j);
                }
            return;
        }

        // Recursively sort halves of dest into src
        int destLow = low;
        int destHigh = high;
        low += off;
        high += off;
        int mid = (low + high) >>> 1;
        mergeSort(dest, src, low, mid, -off, comparator);
        mergeSort(dest, src, mid, high, -off, comparator);

        // If list is already sorted, just copy from src to dest.  This is an
        // optimization that results in faster sorts for nearly ordered lists.
        if (comparator.compare(src[mid - 1], src[mid]) <= 0) {
            System.arraycopy(src, low, dest, destLow, length);
            return;
        }

        // Merge sorted halves (now in src) into dest
        for (int i = destLow, p = low, q = mid; i < destHigh; i++) {
            if (q >= high || p < mid && comparator.compare(src[p], src[q]) <= 0)
                dest[i] = src[p++];
            else
                dest[i] = src[q++];
        }

    }

}
