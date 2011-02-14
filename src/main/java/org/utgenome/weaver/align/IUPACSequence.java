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
// IUPACSequenceReader.java
// Since: 2011/02/10
//
// $URL$ 
// $Author$
//--------------------------------------
package org.utgenome.weaver.align;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.OutputStream;

import org.utgenome.UTGBErrorCode;
import org.utgenome.UTGBException;
import org.utgenome.gwt.utgb.client.bio.IUPAC;
import org.utgenome.weaver.align.SAIS.BaseArray;
import org.xerial.util.FileType;
import org.xerial.util.log.Logger;

/**
 * Reader of the IUPACSequence
 * 
 * @author leo
 * 
 */
public class IUPACSequence implements BaseArray
{
    private static Logger    _logger = Logger.getLogger(IUPACSequence.class);

    private SequenceBoundary binaryInfo;
    private byte[]           seq;
    private int              numBases;

    public IUPACSequence(File iupacFile, int numBases) throws UTGBException {
        initSeq(numBases, iupacFile);
    }

    public IUPACSequence(File iupacFile) throws UTGBException {

        File silkIndexFile = SequenceBoundary.getFileName(FileType.removeFileExt(iupacFile.getPath()));
        binaryInfo = SequenceBoundary.loadSilk(silkIndexFile);

        if (binaryInfo == null)
            throw new UTGBException(UTGBErrorCode.INVALID_INPUT, "failed to load index file");

        initSeq(binaryInfo.totalSize, iupacFile);
    }

    private void initSeq(int numBases, File iupacFile) throws UTGBException {
        this.numBases = numBases;
        int byteSize = (numBases / 2) + (numBases & 0x01);
        this.seq = new byte[byteSize];
        try {

            FileInputStream seqIn = new FileInputStream(iupacFile);
            try {
                _logger.info("loading " + iupacFile);
                int read = seqIn.read(seq, 0, byteSize);
            }
            finally {
                seqIn.close();
            }
        }
        catch (IOException e) {
            throw UTGBException.convert(e);
        }

    }

    public int size() {
        return numBases;
    }

    public IUPAC getIUPAC(int index) {
        byte code = (byte) ((seq[index >> 1] >>> (1 - (index & 1)) * 4) & 0x0F);
        return IUPAC.decode(code);
    }

    public void setIUPAC(int index, IUPAC val) {
        int pos = index / 2;
        int offset = index % 2;
        byte code = (byte) val.bitFlag;

    }

    public IUPAC[] toArray() {
        IUPAC[] result = new IUPAC[size()];
        for (int i = 0; i < result.length; ++i) {
            result[i] = getIUPAC(i);
        }
        return result;
    }

    public void reverse(OutputStream out) throws IOException {
        IUPACSequenceWriter encoder = new IUPACSequenceWriter(out);
        for (int i = binaryInfo.totalSize - 1; i >= 0; --i) {
            encoder.append(this.getIUPAC(i));
        }
        encoder.close();
    }

    @Override
    public int get(int i) {
        return getIUPAC(i).bitFlag;
    }

    @Override
    public void set(int i, int val) {
        throw new UnsupportedOperationException("set");
    }

    @Override
    public int update(int i, int val) {
        throw new UnsupportedOperationException("update");
    }

}
