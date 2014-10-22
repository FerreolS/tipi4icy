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

import icy.gui.frame.progress.AnnounceFrame;
import icy.sequence.Sequence;
import icy.sequence.SequenceEvent;
import icy.sequence.SequenceListener;
import mitiv.deconv.DeconvUtils;
import mitiv.deconv.Deconvolution;
import mitiv.utils.CommonUtils;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.*;
/**
 * EzPlug interface to get the choices of the user
 * 
 * Full CODE see EzPlugTutorial
 * 
 * @author Leger Jonathan
 *
 */
public class MitivDeconvolution extends EzPlug implements EzStoppable,SequenceListener,Block
{

    JSlider slider;

    String wiener = "Wiener";
    String quad = "Quadratic";
    String cg = "CG";

    String normal = "Normal";
    String corrected = "Corrected";
    String colormap = "Colormap";
    String correctColormap = "Corrected+Colormap";

    //Mydata
    EzVarText options = new EzVarText("Regularization", new String[] { wiener,quad,cg}, 0, false);
    EzVarText correction = new EzVarText("Output", new String[] { normal,corrected,colormap,correctColormap}, 0, false);
    EzVarBoolean  varBoolean = new EzVarBoolean("Is PSF splitted ?", false);
    EzVarSequence sequencePSF = new EzVarSequence("PSF");
    EzVarSequence sequenceImage = new EzVarSequence("Image");
    
    EzVarDouble eZcoef = new EzVarDouble("Padding multiplication", 1.0, 10, 0.1);

    EzVarInteger advOne = new EzVarInteger("Complex1");
    EzVarInteger advTwo = new EzVarInteger("Complex2");
    EzVarInteger advThree = new EzVarInteger("Complex3");

    EzVarBoolean advancedOptions = new EzVarBoolean("Show advanced options", false);

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
    }

    public void updateProgressBarMessage(String msg){
        if (!isHeadLess()) {
            getUI().setProgressBarMessage(msg);
        }
    }

    public static double sliderToRegularizationWeight(int slidervalue) {
        return Math.exp(muAlpha + muBeta*slidervalue);
    }

    /*
     * Yes I'm using == to compare 2 strings and yes this is what I want,
     * and yes it's working because getValue() return a string object that I compare
     * with himself
     */
    private int chooseCorrection(){
        if (correction.getValue() == normal) {
            return CommonUtils.SCALE;
        } else if(correction.getValue() == corrected){
            return CommonUtils.SCALE_CORRECTED;
        } else if (correction.getValue() == colormap) {
            return CommonUtils.SCALE_COLORMAP;
        } else if (correction.getValue() == correctColormap){
            return CommonUtils.SCALE_CORRECTED_COLORMAP;
        } else {
            throw new IllegalArgumentException();
        }
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
    private BufferedImage firstJob(int job){
        thread = new ThreadCG(this);
        thread.start();
        boolean isSplitted = varBoolean.getValue();
        switch (job) {
        //First value correspond to next job with alpha = 0, not all are equal to 1
        case DeconvUtils.JOB_WIENER: 
            return (deconvolution.firstDeconvolution(muMin, isSplitted)).get(0);
        case DeconvUtils.JOB_QUAD:
            return (deconvolution.firstDeconvolutionQuad(muMin, isSplitted).get(0));
        case DeconvUtils.JOB_CG:
            return (deconvolution.firstDeconvolutionCG(muMin, isSplitted).get(0));
        default:
            throw new IllegalArgumentException("Invalid Job");
        }
    }

    public BufferedImage nextJob(int slidervalue, int job){
        double mu = sliderToRegularizationWeight(slidervalue);
        if (!isHeadLess()) {
            updateLabel(mu);
        }
        double mult = 1E9; //HACK While the data uniformization is not done...
        switch (job) {
        case DeconvUtils.JOB_WIENER:
            return (deconvolution.nextDeconvolution(mu).get(0));

        case DeconvUtils.JOB_QUAD:
            return (deconvolution.nextDeconvolutionQuad(mu*mult).get(0));

        case DeconvUtils.JOB_CG:
            return (deconvolution.nextDeconvolutionCG(mu*mult).get(0));

        default:
            throw new IllegalArgumentException("Invalid Job");
        }
    }

    private void firstJob3D(int job){
        thread = new ThreadCG(this);
        thread.compute3D();
        thread.start();
        boolean isSplitted = false;
        ArrayList<BufferedImage>tmp;
        switch (job) {
        //First value correspond to next job with alpha = 0, not all are equal to 1
        case DeconvUtils.JOB_WIENER: 
            tmp = deconvolution.firstDeconvolution(muMin,Deconvolution.PROCESSING_3D,isSplitted);
            for (int i = 0; i < tmp.size(); i++) {
                myseq.setImage(0, i, tmp.get(i));
            }
            break;
        case DeconvUtils.JOB_QUAD:
            tmp = deconvolution.firstDeconvolutionQuad(muMin,Deconvolution.PROCESSING_3D,isSplitted);
            for (int i = 0; i < tmp.size(); i++) {
                myseq.setImage(0, i, tmp.get(i));
            }
            break;
        case DeconvUtils.JOB_CG:
            tmp = deconvolution.firstDeconvolutionCG(muMin,Deconvolution.PROCESSING_3D,isSplitted);
            for (int i = 0; i < tmp.size(); i++) {
                myseq.setImage(0, i, tmp.get(i));
            }
            break;
        default:
            throw new IllegalArgumentException("Invalid Job");
        }
    }

    public void nextJob3D(int slidervalue, int job){
        double mu = sliderToRegularizationWeight(slidervalue);
        if (!isHeadLess()) {
            updateLabel(mu);
        }
        double mult = 1E9; //HACK While the data uniformization is not done...
        ArrayList<BufferedImage>tmp;
        switch (job) {
        case DeconvUtils.JOB_WIENER:
            tmp = deconvolution.nextDeconvolution(mu,Deconvolution.PROCESSING_3D);
            for (int i = 0; i < tmp.size(); i++) {
                myseq.setImage(0, i, tmp.get(i));
            }
            break;
        case DeconvUtils.JOB_QUAD:
            tmp = deconvolution.nextDeconvolutionQuad(mu*mult,Deconvolution.PROCESSING_3D);
            for (int i = 0; i < tmp.size(); i++) {
                myseq.setImage(0, i, tmp.get(i));
            }
            break;
        case DeconvUtils.JOB_CG:
            tmp = deconvolution.nextDeconvolutionCG(mu*mult,Deconvolution.PROCESSING_3D);
            for (int i = 0; i < tmp.size(); i++) {
                myseq.setImage(0, i, tmp.get(i));
            }
            break;
        default:
            throw new IllegalArgumentException("Invalid Job");
        }
    }
    
    @Override
    protected void initialize()
    {
        slider = new JSlider(0, 100, 0);
        slider.setEnabled(false);  
        label = new JLabel("                     ");


        addEzComponent(sequencePSF);
        addEzComponent(varBoolean);
        addEzComponent(sequenceImage);
        addEzComponent(options);
        addEzComponent(correction);
        addEzComponent(eZcoef);
        addComponent(slider);
        addComponent(label);
        advancedOptions.addVisibilityTriggerTo(advOne, true);
        advancedOptions.addVisibilityTriggerTo(advTwo, true);
        advancedOptions.addVisibilityTriggerTo(advThree, true);
        addEzComponent(advancedOptions);
        addEzComponent(advOne);
        addEzComponent(advTwo);
        addEzComponent(advThree);

    }

    public void updateImage(BufferedImage buffered, int value){
        if (isHeadLess()) {
            myseq.setImage(0, 0, buffered); 
            output.setValue(myseq);
        } else {
            myseq.setName(options.getValue()+" "+correction.getValue()+" "+value);
            myseq.setImage(0, 0, buffered); 
        }

    }

    @Override
    protected void execute()
    {
        correct = chooseCorrection();
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
            //try {
                Sequence seqIm = sequenceImage.getValue();
                Sequence seqPsf = sequencePSF.getValue();
                if (seqIm.getSizeZ() == 1 && seqPsf.getSizeZ() == 1) {
                    deconvolution = new Deconvolution(seqIm.getFirstNonNullImage(), seqPsf.getFirstNonNullImage(),correct);
                    deconvolution.setPaddingCoefficient(eZcoef.getValue());
                    myseq = new Sequence();
                    myseq.addImage(0,firstJob(job));
                    myseq.addListener(this); 
                    myseq.setName("");
                    if (isHeadLess()) {
                        double value = valueBlock.getValue();
                        updateImage(nextJob((int)value, job), (int)value);
                    } else {
                        addSequence(myseq);
                        slider.setEnabled(true);
                        slider.addChangeListener(new ChangeListener(){
                            public void stateChanged(ChangeEvent event){
                                //getUI().setProgressBarMessage("Computation in progress");
                                int sliderValue =(((JSlider)event.getSource()).getValue());
                                updateProgressBarMessage("Computing");
                                thread.prepareNextJob(sliderValue, job);
                                //OMEXMLMetadataImpl metaData = new OMEXMLMetadataImpl();
                                //myseq.setMetaData(metaData);
                                //updateImage(buffered, tmp);
                            }
                        });  
                        //Beware, need to be called at the END
                        slider.setValue(0);
                    }
                } else {
                    deconvolution = new Deconvolution(seqIm.getAllImage(), seqPsf.getAllImage(),correct);
                    deconvolution.setPaddingCoefficient(eZcoef.getValue());
                    myseq = new Sequence();
                    myseq.addListener(this); 
                    myseq.setName("");
                    firstJob3D(job);
                    if (isHeadLess()) {
                        //TODO pour 3D verifier
                        double value = valueBlock.getValue();
                        updateImage(nextJob((int)value, job), (int)value);
                    } else {
                        addSequence(myseq);
                        slider.setEnabled(true);
                        slider.setValue(0);
                        slider.addChangeListener(new ChangeListener(){
                            public void stateChanged(ChangeEvent event){
                                int sliderValue =(((JSlider)event.getSource()).getValue());
                                updateProgressBarMessage("Computing");
                                thread.prepareNextJob(sliderValue, job);
                            }
                        });
                    }
                }
            //} catch (Exception e) {
            //    new AnnounceFrame("Oops, Error: "+e.getMessage());
            //}
        }
    }

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

    @Override
    public void sequenceClosed(Sequence sequence) {
        if (!isHeadLess()) {
            slider.setEnabled(false);
        }
    }

    public int getOutputValue(){
        return deconvolution.getOuputValue();
    }

    @Override
    public void declareInput(VarList inputMap) {
        inputMap.add(sequencePSF.getVariable());
        inputMap.add(sequenceImage.getVariable());
        inputMap.add(options.getVariable());
        inputMap.add(correction.getVariable());
        inputMap.add(valueBlock.getVariable());
    }

    @Override
    public void declareOutput(VarList outputMap) {
        outputMap.add(output.getVariable());

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