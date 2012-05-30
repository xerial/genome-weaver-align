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
// Strand.java
// Since: 2011/04/21
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.db;

public enum Strand {
	PLUS('+'), MINUS('-');

	public final char symbol;

	private Strand(char symbol) {
		this.symbol = symbol;
	}

	@Override
	public String toString() {
		return Character.toString(symbol);
	}

	public static Strand toStrand(char plusOrMinus) {
		if (plusOrMinus == '+')
			return Strand.PLUS;
		else
			return Strand.MINUS;
	}

	public static Strand toStrand(String plusOrMinus) {
		if (plusOrMinus.equals("+"))
			return Strand.PLUS;
		else
			return Strand.MINUS;
	}

}
