package org.utgenome.weaver.align.strategy;

import java.util.ArrayList;
import java.util.PriorityQueue;

import org.utgenome.gwt.utgb.client.bio.IUPAC;
import org.utgenome.weaver.align.AlignmentSA;
import org.utgenome.weaver.align.AlignmentScoreConfig;
import org.utgenome.weaver.align.BWAlign.FMIndexOnGenome;
import org.utgenome.weaver.align.CharacterCount;
import org.utgenome.weaver.align.FMIndex;
import org.utgenome.weaver.align.IUPACSequence;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.SuffixInterval;
import org.utgenome.weaver.align.record.RawRead;
import org.utgenome.weaver.align.record.ReadSequence;
import org.xerial.util.ObjectHandler;

/**
 * Alignment process using FM-Index.
 * 
 * This class implements BWA-style alignment, breadth-first search with pruning
 * the edges.
 * 
 * @author leo
 * 
 */
public class BWAStrategy
{
    private final FMIndexOnGenome      fmIndex;

    private final int                  numMismatchesAllowed = 1;
    private final AlignmentScoreConfig config               = new AlignmentScoreConfig();
    private final ArrayList<IUPAC>     lettersInGenome      = new ArrayList<IUPAC>();

    public BWAStrategy(FMIndexOnGenome fmIndex) {
        this.fmIndex = fmIndex;

        CharacterCount C = fmIndex.fmIndexF.getCharacterCount();
        for (IUPAC base : IUPAC.values()) {
            if (base == IUPAC.None)
                continue;

            if (C.getCount(base) > 0) {
                lettersInGenome.add(base);
            }
        }
    }

    public void align(RawRead read, ObjectHandler<AlignmentSA> out) throws Exception {

        if (ReadSequence.class.isAssignableFrom(read.getClass())) {
            ReadSequence s = ReadSequence.class.cast(read);
            align(read.name(), new IUPACSequence(s.seq), out);
        }
    }

    /**
     * 
     * @param seq
     * @param cursor
     * @param numMismatchesAllowed
     * @param si
     * @throws Exception
     */
    public void align(String queryName, IUPACSequence seq, ObjectHandler<AlignmentSA> out) throws Exception {

        PriorityQueue<AlignmentSA> alignmentQueue = new PriorityQueue<AlignmentSA>();

        int bestScore = -1;

        int minScore = (int) seq.textSize() * config.matchScore - config.mismatchPenalty * numMismatchesAllowed;
        minScore = Math.max(minScore, (int) (seq.textSize() * 0.5 * config.matchScore));

        alignmentQueue.add(AlignmentSA.initialState(queryName, seq, Strand.FORWARD, fmIndex.fmIndexR.textSize()));
        alignmentQueue.add(AlignmentSA.initialState(queryName, seq.complement(), Strand.REVERSE,
                fmIndex.fmIndexF.textSize()));

        while (!alignmentQueue.isEmpty()) {

            AlignmentSA current = alignmentQueue.poll();
            if (current.numMismatches > numMismatchesAllowed) {
                continue;
            }

            if (current.wordIndex >= seq.textSize()) {
                if (current.alignmentScore >= minScore) {
                    bestScore = current.alignmentScore;
                    out.handle(current);
                }
                continue;
            }

            // Search for deletion
            alignmentQueue.add(current.extendWithDeletion(config));

            IUPAC currentBase = current.common.query.getIUPAC(current.wordIndex);
            // Traverse for each A, C, G, T, ... etc.
            for (IUPAC nextBase : lettersInGenome) {
                FMIndex fm = current.strand == Strand.FORWARD ? fmIndex.fmIndexR : fmIndex.fmIndexF;
                SuffixInterval next = fm.backwardSearch(nextBase, current.suffixInterval);
                if (next.isValidRange()) {
                    // Search for insertion
                    if (current.wordIndex > 0 && current.wordIndex < seq.textSize() - 2) {
                        alignmentQueue.add(current.extendWithInsertion(config));
                    }
                    if ((nextBase.bitFlag & currentBase.bitFlag) != 0) {
                        // match
                        alignmentQueue.add(current.extendWithMatch(next, config));
                    }
                    else {
                        // mismatch
                        alignmentQueue.add(current.extendWithMisMatch(next, config));
                    }
                }
            }
        }

    }
}