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

import java.awt.image.BufferedImage;
import java.text.DecimalFormat;
import java.util.ArrayList;

import javax.swing.JLabel;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import mitiv.array.ShapedArray;
import mitiv.deconv.DeconvUtils;
import mitiv.deconv.Deconvolution;
import mitiv.utils.CommonUtils;
import icy.gui.frame.progress.AnnounceFrame;
import icy.image.IcyBufferedImage;
import icy.sequence.Sequence;
import icy.sequence.SequenceEvent;
import icy.sequence.SequenceListener;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;
import plugins.adufour.ezplug.EzVar;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarListener;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.ezplug.EzVarText;
import plugins.mitiv.io.IcyBufferedImageUtils;

/**
 * EzPlug interface to get the choices of the user
 * 
 * Full CODE see EzPlugTutorial
 * 
 * @author Leger Jonathan
 *
 */
public class MitivWiener extends EzPlug implements EzStoppable, SequenceListener, Block, EzVarListener<String>
{

    JSlider slider;

    String wiener = "Wiener";
    String quad = "Quadratic";
    String cg = "CG";

    //Mydata
    EzVarText options = new EzVarText("Regularization", new String[] { wiener,quad,cg}, 0, false);
    EzVarBoolean  varBoolean = new EzVarBoolean("Is PSF splitted ?", false);
    EzVarSequence sequencePSF = new EzVarSequence("PSF");
    EzVarSequence sequenceImage = new EzVarSequence("Image");
    
    EzVarDouble eZcoef = new EzVarDouble("Padding multiplication", 1.0, 10, 0.1);

    Sequence myseq;
    
    JLabel label;
    int job;
    int correct;
    Deconvolution deconvolution;

    ThreadCG thread;

    //Block
    private EzVarSequence output = new EzVarSequence("Output");
    private EzVarInteger valueBlock = new EzVarInteger("Value");
    //End block
    private final static double muMin = 1e-12;
    private final static double muMax = 1e1;
    private final static double muAlpha = Math.log(muMin);
    private final static double muBeta = Math.log(muMax/muMin)/1e2;

    private void updateLabel(double val){
        DecimalFormat df = new DecimalFormat("#.####");
        label.setText( "Actual Value : "+df.format(val));
        //If in the future we want the output image to have some name
        //myseq.setName(options.getValue()+" "+df.format(val));
    }

    public void updateProgressBarMessage(String msg){
        if (!isHeadLess()) {
            getUI().setProgressBarMessage(msg);
        }
    }

    public static double sliderToRegularizationWeight(int slidervalue) {
        return Math.exp(muAlpha + muBeta*slidervalue);
    }

    private int chooseJob(){
        if (options.getValue().equals(wiener)) {
            return DeconvUtils.JOB_WIENER;
        } else if (options.getValue().equals(quad)) {
            return DeconvUtils.JOB_QUAD;
        } else if (options.getValue().equals(cg)) {
            return DeconvUtils.JOB_CG;
        }else{
            throw new IllegalArgumentException("Invalid Job");
        }
    }

    //FIXME return arraylist bufferedImage and use this for both 3D and 1D
    private IcyBufferedImage firstJob(int job){
        thread = new ThreadCG(this);
        thread.start();
        boolean isSplitted = varBoolean.getValue();
        ShapedArray tmp;
        switch (job) {
        //First value correspond to next job with alpha = 0, not all are equal to 1
        case DeconvUtils.JOB_WIENER:
            tmp = deconvolution.firstDeconvolution(muMin, isSplitted);
            break;
        case DeconvUtils.JOB_QUAD:
            tmp = deconvolution.firstDeconvolutionQuad(muMin, isSplitted);
            break;
        case DeconvUtils.JOB_CG:
            tmp = deconvolution.firstDeconvolutionCG(muMin, isSplitted);
            break;
        default:
            throw new IllegalArgumentException("Invalid Job");
        }
        return IcyBufferedImageUtils.arrayToImage(tmp).get(0);
    }

    public IcyBufferedImage nextJob(int slidervalue, int job){
        double mu = sliderToRegularizationWeight(slidervalue);
        if (!isHeadLess()) {
            updateLabel(mu);
        }
        ShapedArray tmp;
        double mult = 1E9; //HACK While the data uniformization is not done...
        switch (job) {
        case DeconvUtils.JOB_WIENER:
            tmp =deconvolution.nextDeconvolution(mu);
            break;
        case DeconvUtils.JOB_QUAD:
            tmp = deconvolution.nextDeconvolutionQuad(mu*mult);
            break;
        case DeconvUtils.JOB_CG:
            tmp = deconvolution.nextDeconvolutionCG(mu*mult);
            break;
        default:
            throw new IllegalArgumentException("Invalid Job");
        }
        return IcyBufferedImageUtils.arrayToImage(tmp).get(0);
    }

    private void firstJob3D(int job){
        thread = new ThreadCG(this);
        thread.compute3D();
        thread.start();
        boolean isSplitted = false;
        ShapedArray tmp;
        switch (job) {
        //First value correspond to next job with alpha = 0, not all are equal to 1
        case DeconvUtils.JOB_WIENER: 
            tmp = deconvolution.firstDeconvolution(muMin,Deconvolution.PROCESSING_3D,isSplitted);
            break;
        case DeconvUtils.JOB_QUAD:
            tmp = deconvolution.firstDeconvolutionQuad(muMin,Deconvolution.PROCESSING_3D,isSplitted);
            break;
        case DeconvUtils.JOB_CG:
            tmp = deconvolution.firstDeconvolutionCG(muMin,Deconvolution.PROCESSING_3D,isSplitted);

            break;
        default:
            throw new IllegalArgumentException("Invalid Job");
        }
        ArrayList<IcyBufferedImage> list = IcyBufferedImageUtils.arrayToImage(tmp);
        for (int i = 0; i < list.size(); i++) {
            myseq.setImage(0, i, list.get(i));
        }
    }

    public void nextJob3D(int slidervalue, int job){
        double mu = sliderToRegularizationWeight(slidervalue);
        if (!isHeadLess()) {
            updateLabel(mu);
        }
        double mult = 1E10; //HACK While the data uniformization is not done...
        ShapedArray tmp;
        switch (job) {
        case DeconvUtils.JOB_WIENER:
            tmp = deconvolution.nextDeconvolution(mu,Deconvolution.PROCESSING_3D);
            break;
        case DeconvUtils.JOB_QUAD:
            tmp = deconvolution.nextDeconvolutionQuad(mu*mult,Deconvolution.PROCESSING_3D);
            break;
        case DeconvUtils.JOB_CG:
            tmp = deconvolution.nextDeconvolutionCG(mu*mult,Deconvolution.PROCESSING_3D);
            break;
        default:
            throw new IllegalArgumentException("Invalid Job");
        }
        ArrayList<IcyBufferedImage> list = IcyBufferedImageUtils.arrayToImage(tmp);
        for (int i = 0; i < list.size(); i++) {
            myseq.setImage(0, i, list.get(i));
        }
    }
    
    @Override
    protected void initialize()
    {
        slider = new JSlider(0, 100, 0);
        slider.setEnabled(false);  
        label = new JLabel("                     ");

        sequencePSF.setToolTipText(ToolTipText.sequencePSF);
        sequenceImage.setToolTipText(ToolTipText.sequenceImage);
        options.setToolTipText(ToolTipText.textMethod);
        eZcoef.setToolTipText(ToolTipText.doublePadding);
        slider.setToolTipText(ToolTipText.deconvolutionSlider);
        varBoolean.setToolTipText(ToolTipText.booleanPSFSplitted);
        
        addEzComponent(sequencePSF);
        addEzComponent(varBoolean);
        addEzComponent(sequenceImage);
        addEzComponent(options);
        options.addVarChangeListener(this);
        addEzComponent(eZcoef);
        addComponent(slider);
        addComponent(label);

    }

    public void updateImage(BufferedImage buffered, int value){
        if (isHeadLess()) {
            myseq.setImage(0, 0, buffered); 
            output.setValue(myseq);
        } else {
            myseq.setName(options.getValue()+" "+value);
            myseq.setImage(0, 0, buffered); 
        }

    }

    @Override
    protected void execute()
    {
        correct = CommonUtils.SCALE;
        job = chooseJob();
        if (myseq != null) {
            myseq.close();
        }
        //If there is a missing parameter we notify the user with the missing parameter as information
        if(sequenceImage.getValue() == null || sequencePSF.getValue() == null){
            String message = "You have forgotten to give ";
            String messageEnd = "";
            if (sequenceImage.getValue() == null) {
                messageEnd = messageEnd.concat("the image ");
            }
            if(sequencePSF.getValue() == null) {
                if (sequenceImage.getValue() == null) {
                    messageEnd = messageEnd.concat("and ");
                }
                messageEnd = messageEnd.concat("a PSF");
            }
            new AnnounceFrame(message+messageEnd);
        }else{
            try {
                Sequence seqIm = sequenceImage.getValue();
                Sequence seqPsf = sequencePSF.getValue();
                //If there is a 2D image and a 2D psf
                myseq = new Sequence();
                myseq.addListener(this); 
                myseq.setName("");
                if (seqIm.getSizeZ() == 1 && seqPsf.getSizeZ() == 1) {
                    deconvolution = new Deconvolution(seqIm.getFirstNonNullImage(), seqPsf.getFirstNonNullImage(),correct);
                    deconvolution.setPaddingCoefficient(eZcoef.getValue());
                    myseq.addImage(0,firstJob(job));	//For 2D data there is a first job then a next job (faster that way)
                } else if(seqIm.getSizeZ() == seqPsf.getSizeZ()) {
                    ShapedArray imgShapped = IcyBufferedImageUtils.imageToArray(seqIm.getAllImage());
                    ShapedArray psfShapped = IcyBufferedImageUtils.imageToArray(seqPsf.getAllImage());
                    deconvolution = new Deconvolution(imgShapped, psfShapped,correct);
                    deconvolution.setPaddingCoefficient(eZcoef.getValue());
                    firstJob3D(job);
                } else {
                	new AnnounceFrame("The PSF and the image should be of same dimensions");
                	return;
                }
                if (isHeadLess()) {
                    double value = valueBlock.getValue();
                    updateImage(nextJob((int)value, job), (int)value);
                } else {
                    addSequence(myseq);
                    slider.setEnabled(true);
                    slider.addChangeListener(new ChangeListener(){
                        public void stateChanged(ChangeEvent event){
                            int sliderValue =(((JSlider)event.getSource()).getValue());
                            updateProgressBarMessage("Computing");
                            thread.prepareNextJob(sliderValue, job);
                        }
                    });
                    slider.setValue(0);
                }
            } catch (Exception e) {
                new AnnounceFrame("Oops, Error: "+ e.getMessage());
                e.printStackTrace();
            }
        }
    }

    //When the plugin is closed or forced to.
    @Override
    public void clean()
    {
        if(myseq != null){
            myseq.close();
        }
        if(thread != null){
            thread.cancel();
            try {
                thread.join();
            } catch (InterruptedException e) {
                System.err.println("Erreur fin Thread "+e);
            }
        }
    }

    @Override
    public void stopExecution()
    {

    }

    @Override
    public void sequenceChanged(SequenceEvent sequenceEvent) {

    }

    //If the user close a sequence in Icy
    @Override
    public void sequenceClosed(Sequence sequence) {
        if (!isHeadLess()) {
            slider.setEnabled(false);
        }
    }

    //Just a getter 
    public int getOutputValue(){
        return deconvolution.getOuputValue();
    }

    //Inputs for the protocols / blocks
    @Override
    public void declareInput(VarList inputMap) {
        inputMap.add("psf", sequencePSF.getVariable());
        inputMap.add("image", sequenceImage.getVariable());
        inputMap.add("options", options.getVariable());
        inputMap.add("value", valueBlock.getVariable());
    }

    //Outputs for the protocols / blocks
    @Override
    public void declareOutput(VarList outputMap) {
        outputMap.add("output",output.getVariable());

    }

    //Listener to watch the options and enable padding only for CG
    @Override
    public void variableChanged(EzVar<String> source, String newValue) {
        if (newValue == cg) {
            eZcoef.setVisible(true);
        }else {
            eZcoef.setValue(1.0);
            eZcoef.setVisible(false);
        }
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