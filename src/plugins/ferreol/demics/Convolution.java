/*
 * This file is part of TiPi (a Toolkit for Inverse Problems and Imaging)
 * developed by the MitiV project.
 *
 * Copyright (c) 2014 the MiTiV project, http://mitiv.univ-lyon1.fr/
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

package plugins.ferreol.demics;


import static plugins.mitiv.io.Icy2TiPi.sequenceToArray;

import icy.gui.frame.progress.FailedAnnounceFrame;
import icy.plugin.interface_.PluginBundled;
import icy.sequence.Sequence;
import mitiv.array.ShapedArray;
import mitiv.base.Shape;
import mitiv.base.Traits;
import mitiv.linalg.shaped.DoubleShapedVectorSpace;
import mitiv.linalg.shaped.FloatShapedVectorSpace;
import mitiv.linalg.shaped.ShapedVector;
import mitiv.linalg.shaped.ShapedVectorSpace;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarChannel;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.mitiv.io.IcyImager;

/**
 * Icy plugin for convolution
 * @author light
 *
 */
public class Convolution extends EzPlug  implements Block, EzStoppable,PluginBundled {

    protected EzVarSequence psf;
    protected EzVarChannel    psfchannel;        // data channel

    private EzVarBoolean  normalizePSF;

    protected ShapedArray imgArray, psfArray;
    protected Sequence   seqPSF;
    private int vectorSpaceType;
    private ShapedVectorSpace dataSpace, objectSpace;
    private EzVarSequence data;
    private EzVarChannel channel;
    private EzVarSequence outputHeadlessImage;
    private Sequence dataSeq;
    private Shape psfShape;
    private Shape dataShape;



    @Override
    protected void initialize()
    {
        data = new EzVarSequence("Image");
        channel = new EzVarChannel("Data channel:", data.getVariable(), false);
        psf = new EzVarSequence("PSF");
        psfchannel = new EzVarChannel("PSF channel:", psf.getVariable(), false);
        normalizePSF = new EzVarBoolean("Normalize PSF (sum=1):", true);
        addEzComponent(data);
        addEzComponent(channel);
        addEzComponent(psf);
        addEzComponent(psfchannel);
        addEzComponent(normalizePSF);


        if (isHeadLess()) {
            outputHeadlessImage = new EzVarSequence("Output Image");
        }
    }

    @Override
    protected void execute()
    {

        dataSeq = data.getValue();
        seqPSF = psf.getValue();

        if (dataSeq == null)
        {
            throwError("An image should be given");
            return;
        }
        if (seqPSF == null)
        {
            throwError("A psf should be given");
            return;
        }


        imgArray =   sequenceToArray(dataSeq,channel.getValue());
        psfArray =   sequenceToArray(seqPSF,psfchannel.getValue() );
        dataShape = imgArray.getShape();
        psfShape = psfArray.getShape();

        if (imgArray.getType() == Traits.DOUBLE ||
                (psfArray != null && psfArray.getType() == Traits.DOUBLE) ) {
            vectorSpaceType = Traits.DOUBLE;
        } else {
            vectorSpaceType = Traits.FLOAT;
        }

        /* Build vector spaces. */
        if (vectorSpaceType == Traits.FLOAT) {
            dataSpace = new FloatShapedVectorSpace(dataShape);
            objectSpace = new FloatShapedVectorSpace(dataShape);

        } else {
            dataSpace = new DoubleShapedVectorSpace(dataShape);
            objectSpace = new DoubleShapedVectorSpace(dataShape);

        }
        mitiv.conv.Convolution H = mitiv.conv.Convolution.build(dataShape,dataSpace,objectSpace);
        ShapedVector imgVector = dataSpace.create(imgArray);
        ShapedVector resultVector = objectSpace.create();
        H.setPSF(psfArray, null, normalizePSF.getValue());
        H.apply( resultVector ,imgVector);

        Sequence seqY = new Sequence();
        IcyImager.show(resultVector.asShapedArray(),seqY,dataSeq.getName()+"*"+seqPSF.getName(),isHeadLess() );
        if (isHeadLess()) {
            if(outputHeadlessImage==null){
                outputHeadlessImage = new EzVarSequence("Output Image");
            }
            outputHeadlessImage.setValue(seqY);
        }
    }

    @Override
    public void clean() {
        // TODO Auto-generated method stub

    }

    /* (non-Javadoc)
     * @see plugins.adufour.blocks.lang.Block#declareInput(plugins.adufour.blocks.util.VarList)
     */
    @Override
    public void declareInput(VarList inputMap) {
        initialize();
        // TODO Auto-generated method stub
        inputMap.add("image", data.getVariable());
        inputMap.add("image", channel.getVariable());
        inputMap.add("PSF", psf.getVariable());
        inputMap.add("PSF", psfchannel.getVariable());

    }

    /* (non-Javadoc)
     * @see plugins.adufour.blocks.lang.Block#declareOutput(plugins.adufour.blocks.util.VarList)
     */
    @Override
    public void declareOutput(VarList outputMap) {
        // TODO Auto-generated method stub
        outputMap.add("output", outputHeadlessImage.getVariable());

    }

    /* (non-Javadoc)
     * @see icy.plugin.interface_.PluginBundled#getMainPluginClassName()
     */
    @Override
    public String getMainPluginClassName() {
        // TODO Auto-generated method stub
        return "SimpleDEMIC";
    }

    /**
     * print error message
     *
     * @param s
     * the error message
     */
    protected void throwError(final String s) {
        if(isHeadLess()){
            throw new IllegalArgumentException(s);
        } else {
            new FailedAnnounceFrame(s);

        }

    }

}

