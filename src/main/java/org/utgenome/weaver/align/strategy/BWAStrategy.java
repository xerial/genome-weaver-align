package org.utgenome.weaver.align.strategy;

import java.util.ArrayList;
import java.util.PriorityQueue;

import org.utgenome.weaver.align.ACGT;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.AlignmentSA;
import org.utgenome.weaver.align.AlignmentScoreConfig;
import org.utgenome.weaver.align.CharacterCount;
import org.utgenome.weaver.align.FMIndexOnGenome;
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
    private final ArrayList<ACGT>      lettersInGenome      = new ArrayList<ACGT>();

    public BWAStrategy(FMIndexOnGenome fmIndex) {
        this.fmIndex = fmIndex;

        CharacterCount C = fmIndex.forwardIndex.getCharacterCount();
        for (ACGT base : ACGT.values()) {
            if (C.getCount(base) > 0) {
                lettersInGenome.add(base);
            }
        }
    }

    public void align(RawRead read, ObjectHandler<AlignmentSA> out) throws Exception {

        if (ReadSequence.class.isAssignableFrom(read.getClass())) {
            ReadSequence s = ReadSequence.class.cast(read);
            align(s, out);
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
    public void align(ReadSequence read, ObjectHandler<AlignmentSA> out) throws Exception {
        ACGTSequence seq = new ACGTSequence(read.seq);

        int minScore = (int) (seq.textSize() - numMismatchesAllowed) * config.matchScore - config.mismatchPenalty
                * numMismatchesAllowed;
        if (minScore < 0)
            minScore = 1;

        PriorityQueue<AlignmentSA> alignmentQueue = new PriorityQueue<AlignmentSA>();
        int bestScore = -1;

        long N = fmIndex.forwardIndex.textSize();
        alignmentQueue.add(AlignmentSA.initialState(read.name, seq, Strand.FORWARD, N));
        alignmentQueue.add(AlignmentSA.initialState(read.name, seq.complement(), Strand.REVERSE, N));

        while (!alignmentQueue.isEmpty()) {

            AlignmentSA current = alignmentQueue.poll();
            if (current.numMismatches > numMismatchesAllowed) {
                continue;
            }

            if (current.wordIndex >= seq.textSize()) {
                if (current.alignmentScore >= minScore) {
                    if (bestScore < current.alignmentScore)
                        bestScore = current.alignmentScore;
                    out.handle(current);
                }
                continue;
            }

            // search for exact match 
            {
                int remainingMM = numMismatchesAllowed - current.numMismatches;
                if (current.wordIndex == 0 || remainingMM == 0) {
                    SuffixInterval si = alignExact(current.common.query, current.strand);
                    if (si != null) {
                        alignmentQueue.add(current.extendWithMatch(si,
                                (int) (current.common.query.textSize() - current.wordIndex), config));
                        continue;
                    }
                    else {
                        if (remainingMM == 0)
                            continue;
                    }
                }
            }

            // Search for deletion
            alignmentQueue.add(current.extendWithDeletion(config));

            ACGT currentBase = current.common.query.getACGT(current.wordIndex);
            // Traverse for each A, C, G, T, ... etc.
            for (ACGT nextBase : lettersInGenome) {
                SuffixInterval next = fmIndex.forwardSearch(current.strand, nextBase, current.suffixInterval);
                if (next.isValidRange()) {
                    // Search for insertion
                    if (current.wordIndex > 0 && current.wordIndex < seq.textSize() - 1) {
                        alignmentQueue.add(current.extendWithInsertion(config));
                    }
                    if ((nextBase.code & currentBase.code) != 0) {
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

    public SuffixInterval alignExact(ACGTSequence seq, Strand strand) {

        SuffixInterval si = new SuffixInterval(0, fmIndex.forwardIndex.textSize());

        for (int wordIndex = 0; wordIndex < seq.textSize(); ++wordIndex) {
            ACGT currentBase = seq.getACGT(wordIndex);
            SuffixInterval next = fmIndex.forwardSearch(strand, currentBase, si);
            if (!next.isValidRange()) {
                return null;
            }
            si = next;
        }

        return si;
    }

    private void addToQueue(PriorityQueue<AlignmentSA> queue, AlignmentSA newState) {

        queue.add(newState);
    }

}