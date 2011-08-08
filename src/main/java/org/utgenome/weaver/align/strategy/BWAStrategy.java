package org.utgenome.weaver.align.strategy;

import java.util.ArrayList;

import org.utgenome.weaver.align.ACGT;
import org.utgenome.weaver.align.ACGTSequence;
import org.utgenome.weaver.align.AlignmentScoreConfig;
import org.utgenome.weaver.align.CharacterCount;
import org.utgenome.weaver.align.FMIndexOnGenome;
import org.utgenome.weaver.align.Strand;
import org.utgenome.weaver.align.SuffixInterval;
import org.utgenome.weaver.align.record.RawRead;
import org.utgenome.weaver.align.record.ReadSequence;
import org.xerial.util.ObjectHandler;
import org.xerial.util.log.Logger;

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
    private static Logger              _logger              = Logger.getLogger(BWAStrategy.class);

    private final FMIndexOnGenome      fmIndex;
    private final int                  numMismatchesAllowed = 1;
    private final AlignmentScoreConfig config               = new AlignmentScoreConfig();
    private final ArrayList<ACGT>      lettersInGenome      = new ArrayList<ACGT>();
    private ObjectHandler<Alignment>   out;

    public BWAStrategy(FMIndexOnGenome fmIndex) {
        this.fmIndex = fmIndex;

        CharacterCount C = fmIndex.forwardIndex.getCharacterCount();
        for (ACGT base : ACGT.values()) {
            if (C.getCount(base) > 0) {
                lettersInGenome.add(base);
            }
        }
    }

    public void align(RawRead read, ObjectHandler<Alignment> out) throws Exception {
        this.out = out;
        if (ReadSequence.class.isAssignableFrom(read.getClass())) {
            ReadSequence s = ReadSequence.class.cast(read);
            align(s);
        }
    }

    public Alignment exactMatch(Alignment aln) {
        while (!aln.isFinished()) {
            SuffixInterval nextSi = aln.nextSi(fmIndex, aln.nextACGT(), aln.si);
            if (nextSi.isEmpty())
                return null;
            aln = aln.extendWithMatch(config, nextSi);
        }
        return aln;
    }

    public void align(RawRead r) throws Exception {

        // TODO PE mapping
        ReadSequence read = (ReadSequence) r;
        _logger.debug("query: " + read.seq);

        ACGTSequence qF = new ACGTSequence(read.seq);

        if (qF.fastCount(ACGT.N, 0, qF.textSize()) > config.maximumEditDistances) {
            // too many Ns in the query sequence
            return;
        }
        ACGTSequence qC = qF.complement();

        AlignmentQueue queue = new AlignmentQueue(config);
        // Set the initial search states
        queue.add(new Alignment(qF, Strand.FORWARD, SearchDirection.Forward, ExtensionType.MATCH, 0, 0,
                Score.initial(), fmIndex.wholeSARange()));
        queue.add(new Alignment(qC, Strand.REVERSE, SearchDirection.Forward, ExtensionType.MATCH, 0, 0,
                Score.initial(), fmIndex.wholeSARange()));

        // Search iteration
        while (!queue.isEmpty()) {
            Alignment c = queue.poll(); // current 

            if (c.isFinished() && c.score.score >= queue.bestScore) {
                report(c);
                queue.bestScore = c.score.score;
                continue;
            }

            SuffixInterval si = c.suffixInterval();
            if (c.isFinished()) {
                continue;
            }

            int upperBound = c.getUpperBoundOfScore(config);

            if (upperBound < queue.bestScore)
                continue; // no need to proceed

            int remainingDist = config.maximumEditDistances
                    - (c.score.numMismatches + c.score.numGapOpens + c.score.numGapExtend);
            if (remainingDist < 0)
                continue;

            if (_logger.isDebugEnabled())
                _logger.debug("[%3d] %s SI:%s ", queue.queue.size(), c, c.suffixInterval());

            if (remainingDist == 0) {
                if (c.extensionType == ExtensionType.MATCH) {
                    // exact match
                    Alignment a = exactMatch(c);
                    if (a != null)
                        queue.add(a);
                }
                continue;
            }

            // Compute the next SA ranges for A, C, G, T, N
            SuffixInterval[] next = new SuffixInterval[ACGT.exceptN.length]; // for A, C, G, T, N
            {
                int i = 0;
                for (ACGT ch : ACGT.exceptN) {
                    next[i++] = c.nextSi(fmIndex, ch, si);
                }
            }

            // Search for indels
            switch (c.extensionType) {
            case MATCH: { // gap open
                if (c.gapOpenIsAllowed(config) && upperBound - config.gapOpenPenalty > queue.bestScore) {
                    // insertion to reference
                    queue.add(c.startInsertion(config));
                    // deletion from reference
                    for (ACGT ch : ACGT.exceptN) {
                        if (!next[ch.code].isEmpty())
                            queue.add(c.startDeletion(config, next[ch.code]));
                    }
                    break;
                }
            }
            case DELETION: { // gap extension
                if (c.gapExtensionIsAllowed(config) && upperBound - config.gapExtensionPenalty > queue.bestScore) {
                    for (ACGT ch : ACGT.exceptN) {
                        if (!next[ch.code].isEmpty())
                            queue.add(c.extendDeletion(config, next[ch.code]));
                    }
                }
                break;
            }
            case INSERTION: { // gap extension
                if (c.gapExtensionIsAllowed(config) && upperBound - config.gapExtensionPenalty > queue.bestScore) {
                    queue.add(c.extendInsertion(config));
                }
                break;
            }
            }

            // Search for mismatches
            for (ACGT ch : ACGT.exceptN) {
                SuffixInterval nextSi = next[ch.code];
                if (nextSi.isEmpty())
                    continue;

                if (ch == c.nextACGT()) {
                    // match
                    queue.add(c.extendWithMatch(config, nextSi));
                }
                else {
                    if (remainingDist > 0) {
                        // mismatch
                        queue.add(c.extendWithMisMatch(config, nextSi));
                    }
                }
            }

        }
        if (_logger.isDebugEnabled())
            _logger.debug("push count: %,d", queue.pushCount);

    }

    void report(Alignment result) throws Exception {
        out.handle(result);
    }

}