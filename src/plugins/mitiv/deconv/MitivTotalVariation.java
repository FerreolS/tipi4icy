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

import java.util.ArrayList;

import mitiv.array.Double1D;
import mitiv.array.Double2D;
import mitiv.array.Double3D;
import mitiv.array.DoubleArray;
import mitiv.array.ShapedArray;
import mitiv.base.Shape;
import mitiv.invpb.ReconstructionJob;
import mitiv.invpb.ReconstructionViewer;
import mitiv.linalg.WeightGenerator;
import mitiv.utils.FFTUtils;
import icy.gui.frame.progress.AnnounceFrame;
import icy.image.IcyBufferedImage;
import icy.sequence.Sequence;
import icy.sequence.SequenceEvent;
import icy.sequence.SequenceListener;
import commands.TotalVariationDeconvolution;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzGroup;
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

public class MitivTotalVariation extends EzPlug implements Block, EzStoppable, SequenceListener, EzVarListener<String> {
    
    /****************************************************/
    /**                 VIEWER UPDATE                  **/
    /****************************************************/
    public class tvViewer implements ReconstructionViewer{
        @Override
        public void display(ReconstructionJob job) {
            setResult();
        }
    }
    /****************************************************/
    /**                 VARIABLES                      **/
    /****************************************************/
    
    TotalVariationDeconvolution tvDec;

    private Sequence sequence;

    private double mu = 0.1;
    private double epsilon = 0.1;
    private double grtol = 0.1;
    private double gatol = 0.0;
    private int maxIter = 50;
    private boolean goodInput = true;
    private boolean psfSplitted = false;
    private boolean computeNew = true;
    private boolean reuse = true;

    private int width = -1;
    private int height = -1;
    private int sizeZ = -1;
    private double coef = 1.0;
    private Shape shape;

    private EzVarSequence sequencePsf = new EzVarSequence("PSF");
    private EzVarSequence sequenceImg = new EzVarSequence("Image");
    private EzVarSequence output = new EzVarSequence("Output"); //In headLess mode only
    private EzVarBoolean eZpsfSplitted = new EzVarBoolean("Is the psf splitted ?", psfSplitted);
    private EzVarDouble eZmu = new EzVarDouble("Mu", 0, Double.MAX_VALUE, 0.1);
    private EzVarDouble eZepsilon = new EzVarDouble("Epsilon", 0, Double.MAX_VALUE, 1);
    private EzVarDouble eZgrtol = new EzVarDouble("grtol", 0, 1, 0.1);
    private EzVarDouble eZcoef = new EzVarDouble("Padding multiplication", 1.0, 10, 0.1);
    private EzVarInteger eZmaxIter = new EzVarInteger("Max Iterations", -1, Integer.MAX_VALUE, 1);
    private EzVarBoolean eZrestart = new EzVarBoolean("Restart with previous result", false);

    private String weightOption1 = new String("None");
    private String weightOption2 = new String("Personnalized weightMap");
    private String weightOption3 = new String("Variance map");
    private String weightOption4 = new String("Computed Variance");

    private EzVarText options = new EzVarText("Options", new String[] { weightOption1, weightOption2, weightOption3, weightOption4}, 0, false);
    private EzVarSequence weightMap = new EzVarSequence("Weight Map");
    private EzVarSequence varianceMap = new EzVarSequence("Variance Map");
    private EzVarBoolean showPixMap = new EzVarBoolean("Pixel Map ?", false);
    private EzVarSequence deadPixel = new EzVarSequence("Dead pixel map");
    private EzVarDouble alpha = new EzVarDouble("Gain e-/lvl");
    private EzVarDouble beta = new EzVarDouble("Readout noise e-/pxl");
    private EzGroup groupWeighting = new EzGroup("Weighting", options, weightMap, varianceMap, alpha, beta, showPixMap, deadPixel);

    IcyBufferedImage img;
    IcyBufferedImage psf;

    /****************************************************/
    /**                 MESSAGE ICY                    **/
    /****************************************************/
    private void message(String info){
        new AnnounceFrame(info);
        goodInput = false;
    }

    /****************************************************/
    /**                 INITIALIZE ICY                 **/
    /****************************************************/
    @Override
    protected void initialize() {
        //Settings all initials values
        eZmu.setValue(mu);
        eZepsilon.setValue(epsilon);
        eZgrtol.setValue(grtol);
        eZmaxIter.setValue(maxIter);
        eZcoef.setValue(coef);
        options.addVarChangeListener(this);

        //Setting visibility to weights parameters
        alpha.setVisible(false);
        beta.setVisible(false);
        deadPixel.setVisible(false);
        showPixMap.addVisibilityTriggerTo(deadPixel, true);

        //Setting all tooltips
        sequencePsf.setToolTipText(ToolTipText.sequencePSF);
        sequenceImg.setToolTipText(ToolTipText.sequenceImage);
        eZpsfSplitted.setToolTipText(ToolTipText.booleanPSFSplitted);
        eZmu.setToolTipText(ToolTipText.doubleMu);
        eZepsilon.setToolTipText(ToolTipText.doubleEpsilon);
        eZgrtol.setToolTipText(ToolTipText.doubleGrtoll);
        eZmaxIter.setToolTipText(ToolTipText.doubleMaxIter);
        eZcoef.setToolTipText(ToolTipText.doublePadding);
        alpha.setToolTipText(ToolTipText.doubleGain);
        beta.setToolTipText(ToolTipText.doubleNoise);
        eZrestart.setToolTipText(ToolTipText.booleanRestart);
        showPixMap.setToolTipText(ToolTipText.sequencePixel);
        weightMap.setToolTipText(ToolTipText.sequenceWeigth);
        varianceMap.setToolTipText(ToolTipText.sequenceVariance);
        
        //Adding all components
        addEzComponent(sequencePsf);
        addEzComponent(sequenceImg);
        addEzComponent(eZpsfSplitted);
        addEzComponent(eZmu);
        addEzComponent(eZepsilon);
        addEzComponent(eZgrtol);
        addEzComponent(eZmaxIter);
        addEzComponent(eZcoef);
        addEzComponent(eZrestart);

        addEzComponent(groupWeighting);
    }

    /****************************************************/
    /**                  RUN PLUGIN                    **/
    /****************************************************/
    @Override
    protected void execute() {
        //Getting all values
        goodInput = true;
        mu = eZmu.getValue();
        epsilon = eZepsilon.getValue();
        grtol = eZgrtol.getValue();
        maxIter = eZmaxIter.getValue();
        psfSplitted = eZpsfSplitted.getValue();
        coef = eZcoef.getValue();
        if (isHeadLess()) {
            reuse = false;
        }else {
            reuse = eZrestart.getValue();
        }

        //Testing epsilon and grtol
        if (mu < 0) {
            message("Regularization level MU must be strictly positive");
        }
        if (epsilon <= 0) {
            message("Threshold level EPSILON must be strictly positive");
        }
        if (grtol <= 0 || grtol >= 1) {
            message("grtol canno't be lower than 0 or greater than 1");
        }
        if (coef < 1 || coef > 3) {
            message("The Padding can not be lower than 1 or have a value greater than 3");
        }
        if (maxIter < -1)  {
            maxIter = -1;
        }

        //Test if we have the image and the psf ...
        if(sequenceImg.getValue() == null || sequencePsf.getValue() == null){
            //If there is a missing parameter we notify the user with the missing parameter as information
            String message = "You have forgotten to give ";
            String messageEnd = "";
            if (sequenceImg.getValue() == null) {
                messageEnd = messageEnd.concat("the image ");
            }
            if(sequencePsf.getValue() == null) {
                if (sequenceImg.getValue() == null) {
                    messageEnd = messageEnd.concat("and ");
                }
                messageEnd = messageEnd.concat("a PSF");
            }
            message(message+messageEnd);
        }else{
            //And if the sizes are matching
            img = sequenceImg.getValue().getFirstNonNullImage();
            psf = sequencePsf.getValue().getFirstNonNullImage();
            //if the user they the psf is splitted and the psf and image are not of the same size
            if (psfSplitted && (img.getWidth() != psf.getWidth() || img.getHeight() != psf.getHeight())) {
                message("The image and the psf should be of same size");
            }
            //if the user make a mistake between psf and image
            if (psf.getWidth() > img.getWidth() || psf.getHeight() > img.getHeight()) {
                message("The psf canno't be larger than the image");
            }
        }
        //Everything seems good we are ready to launch
        if (goodInput) {
            if (reuse && tvDec != null) {
                //If we restart, we reuse the same data and PSF
                tvDec.setRegularizationWeight(mu);
                tvDec.setRegularizationThreshold(epsilon);
                tvDec.setRelativeTolerance(grtol);
                tvDec.setMaximumIterations(maxIter);
                //Computation HERE
                tvDec.deconvolve(shape);
                //Getting the results
                setResult();
                computeNew = true;
            } else {
                //Launching computation
                tvDec = new TotalVariationDeconvolution();
                tvDec.setRegularizationWeight(mu);
                tvDec.setRegularizationThreshold(epsilon);
                tvDec.setRelativeTolerance(grtol);
                tvDec.setAbsoluteTolerance(gatol);
                tvDec.setMaximumIterations(maxIter);
                tvDec.setViewer(new tvViewer());

                // Read the image and the PSF.
                width = img.getWidth();
                height = img.getHeight();
                sizeZ = sequenceImg.getValue().getSizeZ();

                ArrayList<IcyBufferedImage> listImg = sequenceImg.getValue().getAllImage();
                ArrayList<IcyBufferedImage> listPSf= sequencePsf.getValue().getAllImage();
                DoubleArray imgArray, psfArray;
                double[] weight;
                if (listImg.size() == 1) { //2D
                    double[] image = IcyBufferedImageUtils.icyImage3DToArray1D(listImg, width, height, sizeZ, false);
                    double[] psfTmp = IcyBufferedImageUtils.icyImage3DToArray1D(listPSf, psf.getWidth(), psf.getHeight(), sizeZ, false);
                    weight = createWeight(image);
                    shape = Shape.make(FFTUtils.bestDimension((int)(width*coef)), FFTUtils.bestDimension((int)(height*coef)));
                    imgArray =  Double2D.wrap(image, width, height);
                    psfArray =  Double2D.wrap(psfTmp, psf.getWidth(), psf.getHeight());
                } else { //3D
                    double[] image = IcyBufferedImageUtils.icyImage3DToArray1D(listImg, width, height, sizeZ, false);
                    double[] psfTmp = IcyBufferedImageUtils.icyImage3DToArray1D(listPSf, psf.getWidth(), psf.getHeight(), sizeZ, false);
                    weight = createWeight(image);
                    
                    shape = Shape.make(FFTUtils.bestDimension((int)(width*coef)),
                            FFTUtils.bestDimension((int)(height*coef)),
                            FFTUtils.bestDimension((int)(sizeZ*coef)));
                    imgArray =  Double3D.wrap(image, width, height, sizeZ);
                    psfArray =  Double3D.wrap(psfTmp, psf.getWidth(), psf.getHeight(), sizeZ);
                }

                //BEWARE here we change the value to match the new padded image size
                //addImage(weight, "weights", width, height, sizeZ); //Uncomment this to see weights
                width = FFTUtils.bestDimension((int)(width*coef));
                height = FFTUtils.bestDimension((int)(height*coef));
                sizeZ = FFTUtils.bestDimension((int)(sizeZ*coef));

                tvDec.setWeight(weight);
                tvDec.setData(imgArray);
                tvDec.setPsf(psfArray);

                //Computation HERE
                try{
                tvDec.deconvolve(shape);
                }catch(Exception e){
                    System.err.println(e);
                }
                //Getting the results
                setResult();
                computeNew = true;
            }
            tvDec.setResult(tvDec.getResult());
        }
    }

    /****************************************************/
    /**                  UTILS FUNCTIONS               **/
    /****************************************************/
    //Small utils function that will get the sequence and convert it to ShapedArray
    private ShapedArray weightMapToArray(EzVarSequence seq){
        double[] tmp = IcyBufferedImageUtils.icyImage3DToArray1D(seq.getValue().getAllImage(), width, height, sizeZ, false);
        return Double1D.wrap(tmp, tmp.length);
    }

    private double[] createWeight(double[] data){
        WeightGenerator weightGen = new WeightGenerator();
        String newValue = options.getValue();
        ShapedArray array;
        ShapedArray deadPixMap = null;
        if (showPixMap.getValue() && deadPixel.getValue().getFirstNonNullImage() != null) {
            deadPixMap = weightMapToArray(deadPixel);
        }
        if (newValue == weightOption1) {//No options
            //First case, we create an array of 1, that correspond to the image
            double[] weight = new double[width*height*sizeZ];
            for (int i = 0; i < weight.length; i++) {
                weight[i] = 1;
            }
            weightGen.setWeightMap(Double1D.wrap(weight, weight.length));
        } else if (newValue == weightOption2) {//Weight Map
            //If the user give a weight map we convert it and give it to weightGen
            array = weightMapToArray(weightMap);
            weightGen.setWeightMap(array);
        } else if(newValue == weightOption3) {//Variance map
            //If the user give a varianceMap map we give it to weightGen
            array = weightMapToArray(varianceMap);
            weightGen.setWeightMap(array);//Gain + readOut

        } else if(newValue == weightOption4) {
            //In the case of computed variance: we give gain, readout noise, and the image
            array = Double1D.wrap(data, data.length);
            weightGen.setComputedVariance(array, alpha.getValue(), beta.getValue());

        } else {
            throw new IllegalArgumentException("Incorrect argument for weightmap");
        }
        weightGen.setPixelMap(deadPixMap);
        return weightGen.getWeightMap().toDouble().flatten();//FIXME WRONG but Good for now
    }

    //Debug function Will have to be deleted in the future 
    @SuppressWarnings("unused")
    private void addImage(double[] in, String name, int width, int height, int sizeZ){
        Sequence tmpSeq = new Sequence();
        for (int j = 0; j < sizeZ; j++) {
            double[] temp = new double[width*height];
            for (int i = 0; i < width*height; i++) {
                temp[i] = in[i+j*width*height];
            }
            tmpSeq.setImage(0,j, new IcyBufferedImage(width, height, temp));
        }
        tmpSeq.setName(name);
        addSequence(tmpSeq);
    }

    //The function called by the viewer, here we take care of printing the result in headless mode (protocol) or not
    private void setResult(){
        if (sequence == null || computeNew == true) {
            sequence = new Sequence();
            sequence.addListener(this);
            if (isHeadLess()) {
                output.setValue(sequence);
            }else{
                addSequence(sequence);
            }
            computeNew = false;
        }
        sequence.beginUpdate();
        double[] in = tvDec.getResult().flatten();
        for (int j = 0; j < sizeZ; j++) {
            double[] temp = new double[width*height];
            for (int i = 0; i < width*height; i++) {
                temp[i] = in[i+j*width*height];
            }
            sequence.setImage(0,j, new IcyBufferedImage(width, height, temp));
        }
        sequence.setName("TV mu:"+mu+" Iteration:"+tvDec.getIterations()+" grToll: "+tvDec.getRelativeTolerance());
        sequence.endUpdate();
    }

    //When the plugin is closed we try to close/stop everything
    @Override
    public void clean() {
        if (sequence != null) {
            sequence.close();
        }
        if (tvDec != null) {
            tvDec.stop();
        }
    }

    //Sequence listener: if the sequence is changed
    @Override
    public void sequenceChanged(SequenceEvent sequenceEvent) {
    }

    //Sequence listener: if the sequence is closed
    @Override
    public void sequenceClosed(Sequence seq) {
        //computeNew = true;
    }

    //If the value of ezVarText is changed: meaning we want a particular policy for weightmap
    @Override
    public void variableChanged(EzVar<String> source, String newValue) {
        if (newValue == weightOption1) {
            weightMap.setVisible(false);
            varianceMap.setVisible(false);
            alpha.setVisible(false);
            beta.setVisible(false);
        } else if(newValue == weightOption2) {
            weightMap.setVisible(true);
            varianceMap.setVisible(false);
            alpha.setVisible(false);
            beta.setVisible(false);
        } else if(newValue == weightOption3) {
            weightMap.setVisible(false);
            varianceMap.setVisible(true);
            alpha.setVisible(false);
            beta.setVisible(false);
        } else if(newValue == weightOption4) {
            weightMap.setVisible(false);
            varianceMap.setVisible(false);
            alpha.setVisible(true);
            beta.setVisible(true);
        } else {
            throw new IllegalArgumentException("Incorrect argument for weightmap");
        }
    }

    //If the user call the stop button
    @Override
    public void stopExecution() {
        if (tvDec != null) {
            tvDec.stop();
        }
    }

    //The input variable for the protocol
    @Override
    public void declareInput(VarList inputMap) {
        inputMap.add(sequencePsf.getVariable());
        inputMap.add(sequenceImg.getVariable());
        inputMap.add(eZmu.getVariable());
        inputMap.add(eZepsilon.getVariable());
        inputMap.add(eZgrtol.getVariable());
        inputMap.add(eZmaxIter.getVariable());
        inputMap.add(eZcoef.getVariable());
    }

    //The output variable for the protocol
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