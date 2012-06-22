package org.utgenome.weaver.align.strategy;

import org.utgenome.weaver.align.ACGT;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.AlignmentScoreConfig;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.SiSet;
import org.utgenome.weaver.align.record.AlignmentRecord;
import org.xerial.util.ObjectHandler;

/**
 * Represents an alignment state using FM-index
 * 
 * @author leo
 * 
 */
public class BWAState
{
    public final Cursor        cursor;
    public final ExtensionType extensionType;
    public final Score         score;
    public final SiSet         si;

    protected BWAState(Cursor cursor, ExtensionType extensionType, Score score, SiSet si) {
        this.cursor = cursor;
        this.extensionType = extensionType;
        this.score = score;
        this.si = si;
    }

    protected BWAState(BWAState other) {
        this.cursor = other.cursor;
        this.extensionType = other.extensionType;
        this.score = other.score;
        this.si = other.si;
    }

    public boolean gapOpenIsAllowed(AlignmentScoreConfig config, int posInRead, int m) {
        return score.numGapOpens < config.numGapOpenAllowed && posInRead > config.indelEndSkip
                && (m - posInRead) > config.indelEndSkip;
    }

    public boolean gapExtensionIsAllowed(AlignmentScoreConfig config) {
        return score.numMismatches + score.numGapExtend < config.getMaximumEditDistance(cursor.getFragmentLength());
    }

    public boolean mismatchIsAllowed(AlignmentScoreConfig config) {
        return score.numMismatches < config.getMaximumEditDistance(cursor.getFragmentLength());
    }

    @Override
    public String toString() {
        return String.format("%s %s:%s", cursor.toString(), extensionType.name().charAt(0), score.toString());
    }

    public ACGT nextACGT(ACGTSequence[] q) {
        return cursor.nextACGT(q);
    }

    public int getRemainingBases() {
        return cursor.getRemainingBases();
    }

    public int getUpperBoundOfScore(AlignmentScoreConfig config) {
        return score.score + config.matchScore * getRemainingBases();
    }

    public boolean isFinished() {
        return getRemainingBases() <= 0;
    }

    public SiSet nextSi(FMIndexOnGenome fmIndex, ACGT currentBase) {
        return cursor.nextSi(fmIndex, si, currentBase);
    }

    public BWAState extend(ExtensionType type, int extensionLength, Score newScore, SiSet newSi) {

        Cursor next = cursor;
        while (extensionLength > 0) {
            next = cursor.next();
            --extensionLength;
        }
        return new BWAState(next, type, newScore, newSi);
    }

    public BWAState extendWithMatch(AlignmentScoreConfig config, SiSet newSi) {
        return extend(ExtensionType.MATCH, 1, score.extendWithMatch(config), newSi);
    }

    public BWAState extendWithMisMatch(AlignmentScoreConfig config, SiSet newSi) {
        return extend(ExtensionType.MATCH, 1, score.extendWithMismatch(config), newSi);
    }

    public BWAState startInsertion(AlignmentScoreConfig config) {
        return extend(ExtensionType.INSERTION, 1, score.extendWithGapOpen(config), si);
    }

    public BWAState startDeletion(AlignmentScoreConfig config, SiSet newSi) {
        return extend(ExtensionType.DELETION, 0, score.extendWithGapOpen(config), newSi);
    }

    public BWAState extendInsertion(AlignmentScoreConfig config) {
        return extend(ExtensionType.INSERTION, 1, score.extendWithGapExtend(config), si);
    }

    public BWAState extendDeletion(AlignmentScoreConfig config, SiSet newSi) {
        return extend(ExtensionType.DELETION, 0, score.extendWithGapExtend(config), newSi);
    }

    public void toGenomeCoordinate(FMIndexOnGenome fmIndex, ObjectHandler<AlignmentRecord> reporter) throws Exception {
        final int querySize = cursor.getFragmentLength();

        // TODO report alignment
        //        for (long i = si.lowerBound; i < si.upperBound; ++i) {
        //            PosOnGenome p = fmIndex.toGenomeCoordinate(i, querySize, strand);
        //            if (p != null) {
        //                AlignmentRecord rec = new AlignmentRecord();
        //                rec.chr = p.chr;
        //                rec.start = p.pos;
        //                rec.strand = strand;
        //                rec.score = score.score;
        //                rec.numMismatches = score.numMismatches;
        //                // workaround for Picard tools, which cannot accept base character other than ACGT 
        //                rec.querySeq = strand == Strand.FORWARD ? read.toString() : read.reverse().toString();
        //                // TODO obtain read name
        //                rec.readName = read.toString();
        //                rec.end = p.pos + cursorF;
        //                // TODO output correct CIGAR string
        //                //rec.setCIGAR(result.cigar().toCIGARString());
        //                reporter.handle(rec);
        //            }
        //        }

    }

}