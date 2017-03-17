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


import icy.sequence.Sequence;
import mitiv.array.ShapedArray;
import mitiv.base.Shape;
import mitiv.linalg.shaped.DoubleShapedVector;
import mitiv.linalg.shaped.DoubleShapedVectorSpace;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.mitiv.io.Icy2TiPi;

/**
 * Icy plugin for convolution
 * @author light
 *
 */
public class Convolution extends TiPiPlug {

    private EzVarSequence EzVarSequenceImage;
    private EzVarSequence EzVarSequencePSF;

    @Override
    protected void initialize()
    {
        EzVarSequenceImage = new EzVarSequence("Image");
        EzVarSequencePSF = new EzVarSequence("Load PSF");

        super.addEzComponent(EzVarSequenceImage);
        super.addEzComponent(EzVarSequencePSF);
    }

    @Override
    protected void execute()
    {

        Sequence seqImg = EzVarSequenceImage.getValue();
        Sequence seqPSF = EzVarSequencePSF.getValue();

        int w = seqImg.getSizeX();
        int h = seqImg.getSizeY();
        int d = seqImg.getSizeZ();

        Shape myShape;
        if (d != 1) {
            myShape = new Shape(w, h, d);
        } else {
            myShape = new Shape(w, h);
        }

        ShapedArray imgArray =   Icy2TiPi.sequenceToArray(seqImg);
        ShapedArray psfArray =   Icy2TiPi.sequenceToArray(seqPSF );

        DoubleShapedVectorSpace space = new DoubleShapedVectorSpace(myShape);
        DoubleShapedVector xVector = space.create(imgArray);

        DoubleShapedVector y = space.create();
        mitiv.deconv.Convolution        H = mitiv.deconv.Convolution.build(space);
        H.setPSF(psfArray);
        H.apply( y,xVector);

        Sequence seqY = new Sequence();
        show(y.asShapedArray(),seqY,seqImg.getName()+"*"+seqPSF.getName());

    }

    @Override
    public void clean() {
        // TODO Auto-generated method stub

    }



}

