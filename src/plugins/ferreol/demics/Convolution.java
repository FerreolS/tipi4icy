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


import icy.plugin.interface_.PluginBundled;
import icy.sequence.Sequence;
import mitiv.array.ShapedArray;
import mitiv.base.Shape;
import mitiv.linalg.shaped.DoubleShapedVector;
import mitiv.linalg.shaped.DoubleShapedVectorSpace;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;
import plugins.adufour.ezplug.EzVarChannel;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.ezplug.EzVarText;
import plugins.mitiv.io.Icy2TiPi;
import plugins.mitiv.io.IcyImager;

/**
 * Icy plugin for convolution
 * @author light
 *
 */
public class Convolution extends EzPlug  implements Block, EzStoppable,PluginBundled {

    protected EzVarSequence image;
    protected EzVarChannel    imagechannel;        // data channel
    protected EzVarSequence psf;
    protected EzVarChannel    psfchannel;        // data channel

    protected EzVarText       dataSize;       //
    protected EzVarText       outputSize;     // size of the object after padding
    private EzVarInteger    paddingSizeX,paddingSizeY;
    private int psfSizeX=1,psfSizeY=1,psfSizeZ=1;

    protected EzVarSequence   outputHeadlessImage=null;

    @Override
    protected void initialize()
    {
        image = new EzVarSequence("Image");
        imagechannel = new EzVarChannel("Data channel:", image.getVariable(), false);
        psf = new EzVarSequence("PSF");
        psfchannel = new EzVarChannel("PSF channel:", psf.getVariable(), false);
        addEzComponent(image);
        addEzComponent(imagechannel);
        addEzComponent(psf);
        addEzComponent(psfchannel);
        if (isHeadLess()) {
            outputHeadlessImage = new EzVarSequence("Output Image");
        }
    }

    @Override
    protected void execute()
    {

        Sequence seqImg = image.getValue();
        Sequence seqPSF = psf.getValue();

        int w = seqImg.getSizeX();
        int h = seqImg.getSizeY();
        int d = seqImg.getSizeZ();

        Shape myShape;
        if (d != 1) {
            myShape = new Shape(w, h, d);
        } else {
            myShape = new Shape(w, h);
        }

        ShapedArray imgArray =   Icy2TiPi.sequenceToArray(seqImg,imagechannel.getValue());
        ShapedArray psfArray =   Icy2TiPi.sequenceToArray(seqPSF,psfchannel.getValue() );

        DoubleShapedVectorSpace space = new DoubleShapedVectorSpace(myShape);
        DoubleShapedVector xVector = space.create(imgArray);

        DoubleShapedVector y = space.create();
        mitiv.deconv.Convolution        H = mitiv.deconv.Convolution.build(space);
        H.setPSF(psfArray);
        H.apply( y,xVector);

        Sequence seqY = new Sequence();
        IcyImager.show(y.asShapedArray(),seqY,seqImg.getName()+"*"+seqPSF.getName(),isHeadLess() );
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
        inputMap.add("image", image.getVariable());
        inputMap.add("image", imagechannel.getVariable());
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



}

