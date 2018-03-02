/*  Copyright (C) 2017  Ferreol Soulez ferreol.soulez@univ-lyon1.fr
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  */

package plugins.ferreol.demics;

import icy.gui.frame.progress.AnnounceFrame;
import icy.image.IcyBufferedImage;
import icy.plugin.interface_.PluginBundled;
import icy.sequence.Sequence;
import mitiv.array.ShapedArray;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.mitiv.io.Icy2TiPi;

/**
 * This plugin is normalizing all the data so that the sum of all pixels
 * is equal to one.
 * @deprecated
 * @author light
 *
 */
@Deprecated
public class MitivNormalization extends EzPlug implements  EzStoppable, Block, PluginBundled  {

    EzVarSequence image = new EzVarSequence("Image to normalize");
    //Block
    EzVarSequence imageOut = new EzVarSequence("Image Normalized");
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
            ShapedArray tmp1 =  Icy2TiPi.sequenceToArray(seq); //By default we take the first canal
            double[] out = tmp1.toDouble().flatten();
            for (int i = 0; i < out.length; i++) {
                count += out[i];
            }
            //Normalization
            System.out.println("SUM: "+count);
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
            if (isHeadLess()) {
                imageOut.setValue(seqOut);
            } else {
                addSequence(seqOut);
            }
        }
    }

    @Override
    public void clean() {
    }

    @Override
    public void declareInput(VarList inputMap) {
        inputMap.add("imageIn", image.getVariable());
    }

    @Override
    public void declareOutput(VarList outputMap) {
        outputMap.add("imageOut", imageOut.getVariable());
    }

    @Override
    public String getMainPluginClassName() {
        return  SimpleDEMIC.class.getName();
    }



}