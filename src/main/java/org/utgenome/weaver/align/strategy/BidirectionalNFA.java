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
// BidirectionalNFA.java
// Since: Aug 5, 2011
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align.strategy;

import java.util.ArrayList;
import java.util.List;

import org.utgenome.weaver.align.ACGT;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.SuffixInterval;
import org.utgenome.weaver.parallel.Reporter;
import org.xerial.util.Optional;

/**
 * NFA for forward and backward traversal of suffix arrays
 * 
 * @author leo
 * 
 */
public class BidirectionalNFA {

	private final FMIndexOnGenome fmIndex;
	private final ACGTSequence query;
	private final Strand strand;
	private final int m; // query length

	public BidirectionalNFA(FMIndexOnGenome fmIndex, ACGTSequence query,
			Strand strand) {
		this.fmIndex = fmIndex;
		this.query = query;
		this.strand = strand;
		this.m = (int) query.textSize();
	}

	private static class CursorContainer {
		private List<Optional<List<Cursor>>> container;

		public CursorContainer(int m) {
			container = new ArrayList<Optional<List<Cursor>>>(m);
			for (int i = 0; i < m; ++i) {
				container.add(new Optional<List<Cursor>>());
			}
		}

		public boolean hasEntry(int index) {
			return container.get(index).isDefined();
		}

		public List<Cursor> get(int index) {
			return container.get(index).get();
		}

		public void add(int index, Cursor c) {
			Optional<List<Cursor>> l = container.get(index);
			if (l.isUndefined()) {
				l.set(new ArrayList<Cursor>());
			}
			l.get().add(c);
		}

	}

	private static class Cursor {
		public final SuffixInterval si;
		public final int cursorF;
		public final int cursorB;

		public Cursor(SuffixInterval si, int cursorF, int cursorB) {
			this.si = si;
			this.cursorF = cursorF;
			this.cursorB = cursorB;
		}
	}

	public void align(Reporter out) {

		// Row-wise simulation of NFA
		CursorContainer stateF = new CursorContainer(m + 1);
		CursorContainer stateB = new CursorContainer(m + 1);

		// Simulate k=0 (no mismatch)
		// Forward search
		{
			SuffixInterval si = fmIndex.wholeSARange();
			int i = 0, breakPoint = 0;
			stateF.add(0, new Cursor(si, i, 0));
			while (i < m) {
				ACGT ch = query.getACGT(i);
				SuffixInterval nextSi = fmIndex.forwardSearch(strand, ch, si);
				if (!nextSi.isEmpty()) {
					++i;
					stateF.add(i, new Cursor(nextSi, i, breakPoint));
				} else {
					breakPoint = i;
					si = fmIndex.wholeSARange();
					if (breakPoint == 0) {
						// exact match
						// out.emit(result);
					}
				}
			}
		}

		// Backward search from the last match
		if (stateF.hasEntry(m)) {
			for (Cursor each : stateF.get(m)) {
				// compute reverse SA range
				SuffixInterval rSi = fmIndex.wholeSARange();
				for (int i = each.cursorB; i < each.cursorF; ++i) {
					ACGT ch = query.getACGT(i);
					//fmIndex.bidirectionalSearch(strand, nextBase, forwardSi, backwardSi)
					
				}
			}
		}

	}

}
