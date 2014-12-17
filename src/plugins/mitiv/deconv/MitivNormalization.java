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

package plugins.mitiv.deconv;

import icy.gui.frame.progress.AnnounceFrame;
import icy.image.IcyBufferedImage;
import icy.sequence.Sequence;
import mitiv.array.ShapedArray;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.mitiv.io.IcyBufferedImageUtils;

public class MitivNormalization extends EzPlug {

    EzVarSequence image = new EzVarSequence("Image to normalize");
    boolean goodInput = true;
    
    @Override
    protected void initialize() {
        addEzComponent(image);
    }

    private void message(String info){
        new AnnounceFrame(info);
        goodInput = false;
    }
    
    @Override
    protected void execute() {
        Sequence seq = image.getValue();
        if (seq == null) {
            message("You should give a image");
        }
        if (goodInput) {
            IcyBufferedImage img = seq.getFirstNonNullImage();
            int width = img.getWidth();
            int height = img.getHeight();
            int sizeZ = seq.getSizeZ();
            double count = 0.0;
            //Icy to double
            ShapedArray tmp1 =  IcyBufferedImageUtils.imageToArray(seq, 0); //By default we take the first canal
            double[] out = tmp1.toDouble().flatten();
            for (int i = 0; i < out.length; i++) {
                count += out[i];
            }
            //Normalization
            System.out.println(count);
            //Double to icy
            Sequence seqOut = new Sequence();
            for (int j = 0; j < sizeZ; j++) {
                double[] tmp = new double[width*height];
                for (int i = 0; i < tmp.length; i++) {
                    tmp[i] = out[i+j*tmp.length]/count;
                }
                seqOut.setImage(0,j,new IcyBufferedImage(width, height, tmp));
            }
            //Add the sequence to icy
            seqOut.setName("Normalized_"+seq.getName());
            addSequence(seqOut);
        }
    }
    
    @Override
    public void clean() {
    }



}


/*
 * Local Variables:
 * mode: Java
 * tab-width: 8
 * indent-tabs-mode: nil
 * c-basic-offset: 4
 * fill-column: 78
 * coding: utf-8
 * ispell-local-dictionary: "american"
 * End:
 */