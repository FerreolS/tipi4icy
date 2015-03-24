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


import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTextPane;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import mitiv.linalg.shaped.DoubleShapedVector;
import mitiv.utils.MathUtils;
import icy.sequence.Sequence;
import icy.sequence.SequenceEvent;
import icy.sequence.SequenceListener;
import icy.gui.frame.GenericFrame;
import icy.image.IcyBufferedImage;
import icy.type.collection.array.Array1DUtil;
import plugins.adufour.ezplug.*;

/**
 * A demo plugin used to create blurred for a presentation.
 * 
 * Note: This plugin was not updated and may (will) not works.
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
        
        EzButton detailsButton = new EzButton("Help", new ActionListener()
        {
            @Override
            public void actionPerformed(ActionEvent e)
            {
                onDetailsClicked();
            }
        });
        addEzComponent(detailsButton);
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
    
    /**
     * Action performed when the user clicks on the details button
     */
    private void onDetailsClicked()
    {
        String    title   = "Wild Field Fluorescent Microscopy 3D Blind Deconvolution";
        JTextPane message = new JTextPane();
        message.setEditable(false);
        message.setContentType("text/html");
        message.setText(
                "<p>" +
                        "In what follows, <tt>A</tt> represents the input clean sequence, " +
                        "while <tt>B</tt> is the output noisy sequence." +
                        "</p>" +
                        "<h2>Data</h2>" +
                        "<p>" +
                        "1 parameter: <tt>sigma &gt;= 0</tt> (standard deviation of the Gaussian random variables)." +
                        "</p>" +
                        "<p>" +
                        "Exemple deltaX = ni/lambda/4 = 700185 and deltaY = ni/lambda = 2800740" +
                        "zdepth = nz*dz = 10.24 um" +
                        "This noise model enforces <tt>B = A + n</tt>, where <tt>n</tt> is a " +
                        "random sequence such that the samples <tt>n(x,y,z,t,c)</tt> are random " +
                        "independant variables following a Gaussian probability distribution " +
                        "of mean <tt>0</tt> and variance <tt>sigma^2</tt>." +
                        "</p>" +

            "<h2>Microscope and point spread function settings</h2>" +
            "<p>" +
            "No parameter." +
            "</p>" +
            "<p>" +
            "In this model, each output sample <tt>B(x,y,z,t,c)</tt> is generated from a Poisson random " +
            "distribution of intensity <tt>A(x,y,z,t,c)</tt>, which is supposed to be <tt>&gt;=0</tt>. " +
            "If <tt>A(x,y,z,t,c) &lt; 0</tt>, then <tt>B(x,y,z,t,c)</tt> is set to <tt>NaN</tt>." +
            "</p>" +

            "<h2>Weight</h2>" +
            "<p>" +
            "1 parameter: <tt>sigma &gt;= 0</tt> (standard deviation of the Gaussian random variables)." +
            "</p>" +

            "<h2>Optimisation</h2>" +
            "<p>" +
            "1 parameter: <tt>sigma &gt;= 0</tt> (standard deviation of the Gaussian random variables)." +
            "</p>" +

            "<h2>Hints</h2>" +
            "<p>" +
            "<ul>" +
            "<li>Use the Icy colormap in the lookup table at the left to change the colormap model</li>" +
            "<li>Use 3D OrthoViewer to see YX and XZ views. To activate this view mode," +
            "select the Orthoviewer logo in the viewer window (in the dropdown list that defaults to 2D mode).</li>" +
            "</ul>" +
            "</p>"
                );
        Dimension dim = message.getPreferredSize();
        dim.setSize(600, dim.getHeight()+100);
        message.setPreferredSize(dim);
        JScrollPane scroll = new JScrollPane(message);
        dim = scroll.getPreferredSize();
        dim.setSize(600, 500);
        scroll.setPreferredSize(dim);
        GenericFrame infoFrame = new GenericFrame(title, scroll);
        infoFrame.addToDesktopPane();
        infoFrame.setVisible(true);
        infoFrame.requestFocus();
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