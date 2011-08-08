package org.utgenome.weaver.align.strategy;

import org.utgenome.weaver.align.ACGT;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.AlignmentScoreConfig;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.SuffixInterval;

/**
 * Represents an alignment state
 * 
 * @author leo
 * 
 */
public class Alignment
{
    // minimum bit length for this state: 1,1,2,2,16,16,32+8+8+8,64+64 = 222 bit < 28 bytes
    public final ACGTSequence    read;
    public final Strand          strand;       // forward or reverse (complement of the query sequence)
    public final SearchDirection orientation;  // search direction
    public final ExtensionType   extensionType;
    public final int             cursorF;
    public final int             cursorB;
    public final Score           score;
    public final SuffixInterval  si;

    protected Alignment(ACGTSequence read, Strand strand, SearchDirection orientation, ExtensionType extensionType,
            int cursorF, int cursorB, Score score, SuffixInterval si) {
        this.read = read;
        this.strand = strand;
        this.orientation = orientation;
        this.extensionType = extensionType;
        this.cursorF = cursorF;
        this.cursorB = cursorB;
        this.score = score;
        this.si = si;
    }

    protected Alignment(Alignment other) {
        this.read = other.read;
        this.strand = other.strand;
        this.orientation = other.orientation;
        this.extensionType = other.extensionType;
        this.cursorF = other.cursorF;
        this.cursorB = other.cursorB;
        this.score = other.score;
        this.si = other.si;
    }

    public boolean gapOpenIsAllowed(AlignmentScoreConfig config) {
        return score.numGapOpens < config.numGapOpenAllowed;
    }

    public boolean gapExtensionIsAllowed(AlignmentScoreConfig config) {
        return score.numMismatches + score.numGapExtend < config.maximumEditDistances;
    }

    public boolean mismatchIsAllowed(AlignmentScoreConfig config) {
        return score.numMismatches < config.maximumEditDistances;
    }

    @Override
    public String toString() {
        return String.format("%s%s%s(%2d,%2d):%s", strand.symbol, orientation.symbol, extensionType.name().charAt(0),
                cursorF, cursorB, score.toString());
    }

    public ACGT nextACGT() {
        int cursor = orientation.isForward ? cursorF : cursorB - 1;
        return read.getACGT(cursor);
    }

    public int getRemainingBases() {
        int M = (int) read.textSize();
        return (M - cursorF) + cursorB;
    }

    public SuffixInterval suffixInterval() {
        return si;
    }

    public int getUpperBoundOfScore(AlignmentScoreConfig config) {
        return score.score + config.matchScore * getRemainingBases();
    }

    public boolean isFinished() {
        return getRemainingBases() <= 0;
    }

    public SuffixInterval nextSi(FMIndexOnGenome fmIndex, ACGT ch, SuffixInterval si) {
        if (orientation.isForward)
            return fmIndex.forwardSearch(strand, ch, si);
        else
            return fmIndex.backwardSearch(strand, ch, si);
    }

    public Alignment extend(ExtensionType type, int extensionLength, Score newScore, SuffixInterval newSi) {
        int M = (int) read.textSize();
        int nextF = cursorF;
        int nextB = cursorB;
        if (orientation.isForward)
            nextF += extensionLength;
        else
            nextB -= extensionLength;

        Alignment next = null;
        switch (orientation) {
        case Forward:
        case Backward:
            next = new Alignment(read, strand, orientation, type, nextF, nextB, newScore, newSi);
            break;
        case BidirectionalForward:
            if (nextF >= M) {
                // switch to backward search

            }
            break;
        default:
            throw new IllegalStateException("cannot reach here");
        }

        return next;
    }

    public Alignment extendWithMatch(AlignmentScoreConfig config, SuffixInterval newSi) {
        return extend(ExtensionType.MATCH, 1, score.extendWithMatch(config), newSi);
    }

    public Alignment extendWithMisMatch(AlignmentScoreConfig config, SuffixInterval newSi) {
        return extend(ExtensionType.MATCH, 1, score.extendWithMismatch(config), newSi);
    }

    public Alignment startInsertion(AlignmentScoreConfig config) {
        return extend(ExtensionType.INSERTION, 1, score.extendWithGapOpen(config), si);
    }

    public Alignment startDeletion(AlignmentScoreConfig config, SuffixInterval newSi) {
        return extend(ExtensionType.DELETION, 0, score.extendWithGapOpen(config), newSi);
    }

    public Alignment extendInsertion(AlignmentScoreConfig config) {
        return extend(ExtensionType.INSERTION, 1, score.extendWithGapExtend(config), si);
    }

    public Alignment extendDeletion(AlignmentScoreConfig config, SuffixInterval newSi) {
        return extend(ExtensionType.DELETION, 0, score.extendWithGapExtend(config), newSi);
    }

}