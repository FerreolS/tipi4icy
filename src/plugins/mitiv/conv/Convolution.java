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

package plugins.mitiv.conv;


import javax.swing.JLabel;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import mitiv.linalg.shaped.DoubleShapedVector;
import mitiv.utils.MathUtils;
import icy.sequence.Sequence;
import icy.sequence.SequenceEvent;
import icy.sequence.SequenceListener;
import icy.image.IcyBufferedImage;
import icy.type.collection.array.Array1DUtil;
import plugins.adufour.ezplug.*;

/**
 * EzPlug interface to get the choices of the user
 * 
 * Full CODE see EzPlugTutorial
 * 
 * @author Ludovic Beaubras
 *
 */
public class Convolution extends EzPlug implements EzStoppable,SequenceListener, EzVarListener<String>
{
    //Mydata
    EzVarText   varText;
    EzVarText  options;
    EzVarText  kernel;
    EzVarText  noise;
    //EzVarBoolean  varBoolean;
    EzVarFile   varFilePSF;
    EzVarFile   varFileIMAGE;
    EzVarSequence EzVarSequencePSF;
    EzVarSequence EzVarSequenceImage;
    JSlider slider;
    int sliderValue = 3;

    String[] filters = {"no kernel", "average", "disk", "sobel", "prewitt", "kirsh"};
    String noNoise = "no noise";
    String gaussian = "gaussian";

    Sequence myseq;
    Sequence myseqData;
    JLabel label;


    public void updateProgressBarMessage(String msg){
        getUI().setProgressBarMessage(msg);
    }

    public static final int GAUSSIAN = 2;
    public static final int AVERAGE = 3;
    public static final int PREWITT = 4;
    public static final int SOBEL = 5;
    public static final int KIRSH = 6;
    public static final int DISK = 7;
    

    @Override
    protected void initialize()
    {
        EzVarSequenceImage = new EzVarSequence("Image");
        EzVarSequencePSF = new EzVarSequence("Load PSF");
        noise = new EzVarText("Noise", new String[] {noNoise, gaussian}, 0, false);
        slider = new JSlider(0, 100, sliderValue);
        slider.setEnabled(false);  
        label = new JLabel("Value : " + sliderValue);
        slider.addChangeListener(new ChangeListener(){
            public void stateChanged(ChangeEvent event){
                sliderValue =(((JSlider)event.getSource()).getValue());
                label.setText("Value : " + sliderValue);
            }
        });

        super.addEzComponent(EzVarSequenceImage);
        super.addEzComponent(EzVarSequencePSF);
        super.addEzComponent(noise);
        super.addComponent(slider);
        super.addComponent(label);
    }

    @Override
    protected void execute()
    {
        if (myseq != null) {
            myseq.close();
        }
        
        Sequence seqImg = EzVarSequenceImage.getValue();
        Sequence seqPSF = EzVarSequencePSF.getValue();

        int w = seqImg.getSizeX();
        int h = seqImg.getSizeY();
        int d = seqImg.getSizeZ();

        double[] x = new double[w*h*d];
        double[] psf = new double[w*h*d];
        for(int k = 0; k < d; k++)
        {
            Array1DUtil.arrayToDoubleArray(seqImg.getDataXY(0, k, 0), 0, x, k*w*h, w*h, seqImg.isSignedDataType());
            Array1DUtil.arrayToDoubleArray(seqPSF.getDataXY(0, k, 0), 0, psf, k*w*h, w*h, seqImg.isSignedDataType());
        }
        
            double[] psf_shift= MathUtils.fftShift3D(psf, w, h, d);
            DoubleShapedVector y = MathUtils.convolution(x, psf_shift, w, h, d);
            
            Sequence seqY = new Sequence();
            seqY.setName("Y");
            for(int k = 0; k < d; k++)
            {
                seqY.addImage(new IcyBufferedImage(w, h, MathUtils.getArray(y.getData(), w, h, k)));
            }
            addSequence(seqY);

    }

    @Override
    public void clean()
    {
        if(myseq != null){
            myseq.close();
        }
    }

    @Override
    public void stopExecution()
    {

    }

    @Override
    public void sequenceChanged(SequenceEvent sequenceEvent) {

    }

    @Override
    public void sequenceClosed(Sequence sequence) {
        //slider.setEnabled(false);   
    }
    @Override
    public void variableChanged(EzVar<String> source, String newValue) {
        boolean chosen = newValue.compareTo("average") == 0 || newValue.compareTo("disk") == 0;
        slider.setEnabled(chosen);
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