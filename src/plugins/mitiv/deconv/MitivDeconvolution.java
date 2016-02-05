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

import loci.formats.ome.OMEXMLMetadataImpl;
import mitiv.array.Double1D;
import mitiv.array.DoubleArray;
import mitiv.array.ShapedArray;
import mitiv.base.Shape;
import mitiv.base.mapping.DoubleFunction;
import mitiv.invpb.ReconstructionJob;
import mitiv.invpb.ReconstructionViewer;
import mitiv.linalg.WeightGenerator;
import mitiv.utils.FFTUtils;
import mitiv.utils.reconstruction.ReconstructionThread;
import mitiv.utils.reconstruction.ReconstructionThreadToken;
import icy.gui.frame.progress.AnnounceFrame;
import icy.image.IcyBufferedImage;
import icy.sequence.Sequence;
import icy.sequence.SequenceEvent;
import icy.sequence.SequenceListener;
import icy.util.OMEUtil;
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
import plugins.mitiv.reconstruction.TotalVariationJobForIcy;
/**
 * MitivTotalVariation is intended to offer a powerful deconvolution tool with 
 * a nice interface and options for advanced users.
 * 
 * @author light
 *
 */
public class MitivDeconvolution extends EzPlug implements Block, EzStoppable, SequenceListener, EzVarListener<String> {

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

    TotalVariationJobForIcy tvDec;

    private Sequence sequence;

    private double mu = 0.1;
    private double epsilon = 0.1;
    private double grtol = 0.0;
    private double gatol = 0.0;
    private double lowerBound = Double.NEGATIVE_INFINITY;
    private int maxIter = 50;
    private boolean goodInput = true;
    //private boolean psfSplitted = false;
    private boolean computeNew = true;
    private boolean reuse = false;

    private int width = -1;
    private int height = -1;
    private int sizeZ = -1;
    private double coef = 1.0;
    private Shape shape;

    private ReconstructionThreadToken token;
    ReconstructionThread thread;

    private EzVarSequence sequencePsf;
    private EzVarSequence sequenceImg;
    private EzVarSequence lastResult;
    private EzVarSequence output;
    //private EzVarBoolean eZpsfSplitted = new EzVarBoolean("Is the psf splitted?", psfSplitted);
    private EzVarDouble eZmu;
    private EzVarDouble eZepsilon;
    private EzVarInteger eZcoef;
    private EzVarInteger eZmaxIter;
    private EzVarBoolean eZpositivity;

    private String weightOption1;
    private String weightOption2;
    private String weightOption3;
    private String weightOption4;

    private EzVarText options;
    private EzVarSequence weightMap;
    private EzVarSequence varianceMap;
    private EzVarBoolean showPixMap;
    private EzVarSequence deadPixel;
    private EzVarDouble alpha;
    private EzVarDouble beta;
    private EzGroup groupRegularization;
    private EzGroup groupWeighting;
    private EzGroup groupConvergence;

    IcyBufferedImage img;
    IcyBufferedImage psf;
    IcyBufferedImage result;

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

        sequencePsf = new EzVarSequence("Input PSF");
        sequenceImg = new EzVarSequence("Input Image");
        lastResult = new EzVarSequence("Previous result");
        output = new EzVarSequence("Output"); //In headLess mode only
        //private EzVarBoolean eZpsfSplitted = new EzVarBoolean("Is the psf splitted?", psfSplitted);
        eZmu = new EzVarDouble("Regularization level", 0, Double.MAX_VALUE, 0.1);
        eZepsilon = new EzVarDouble("Threshold level", 0, Double.MAX_VALUE, 1);
        eZcoef = new EzVarInteger("Number of lines to add (padding)", 0, 10000, 1);
        eZmaxIter = new EzVarInteger("Max Iterations", -1, Integer.MAX_VALUE, 1);
        eZpositivity = new EzVarBoolean("Enforce nonnegativity", false);

        weightOption1 = new String("None");
        weightOption2 = new String("Personnalized weightMap");
        weightOption3 = new String("Variance map");
        weightOption4 = new String("Computed Variance");

        options = new EzVarText("Options", new String[] { weightOption1, weightOption2, weightOption3, weightOption4}, 0, false);
        weightMap = new EzVarSequence("Weight Map");
        varianceMap = new EzVarSequence("Variance Map");
        showPixMap = new EzVarBoolean("Data Map?", false);
        deadPixel = new EzVarSequence("Dead pixel map");
        alpha = new EzVarDouble("Gain e-/lvl");
        beta = new EzVarDouble("Readout noise e-/pxl");
        groupRegularization = new EzGroup("Regularization", eZmu, eZepsilon);
        groupWeighting = new EzGroup("Weighting", options, weightMap, varianceMap, alpha, beta, showPixMap, deadPixel);
        groupConvergence = new EzGroup("Convergence settings", eZmaxIter);

        lastResult.setNoSequenceSelection();

        //Settings all initials values
        eZmu.setValue(mu);
        eZepsilon.setValue(epsilon);
        eZmaxIter.setValue(maxIter);
        eZcoef.setValue(0);
        options.addVarChangeListener(this);

        //Setting visibility to weights parameters
        alpha.setVisible(false);
        beta.setVisible(false);
        deadPixel.setVisible(false);
        showPixMap.addVisibilityTriggerTo(deadPixel, true);

        //Setting all tooltips
        sequencePsf.setToolTipText(ToolTipText.sequencePSF);
        sequenceImg.setToolTipText(ToolTipText.sequenceImage);
        lastResult.setToolTipText(ToolTipText.booleanRestart);
        //eZpsfSplitted.setToolTipText(ToolTipText.booleanPSFSplitted);
        eZmu.setToolTipText(ToolTipText.doubleMu);
        eZepsilon.setToolTipText(ToolTipText.doubleEpsilon);
        eZmaxIter.setToolTipText(ToolTipText.doubleMaxIter);
        eZcoef.setToolTipText(ToolTipText.doublePadding);
        alpha.setToolTipText(ToolTipText.doubleGain);
        beta.setToolTipText(ToolTipText.doubleNoise);
        eZpositivity.setToolTipText(ToolTipText.booleanPositivity);
        showPixMap.setToolTipText(ToolTipText.sequencePixel);
        weightMap.setToolTipText(ToolTipText.sequenceWeigth);
        varianceMap.setToolTipText(ToolTipText.sequenceVariance);

        //Adding all components
        addEzComponent(sequenceImg);
        addEzComponent(sequencePsf);
        addEzComponent(lastResult);
        //addEzComponent(eZpsfSplitted);
        addEzComponent(groupRegularization);
        addEzComponent(groupConvergence);
        addEzComponent(groupWeighting);
        addEzComponent(eZcoef);
        addEzComponent(eZpositivity);

        token = new ReconstructionThreadToken(new double[]{mu,epsilon,gatol,grtol});
        thread = new ReconstructionThread(token);
        thread.start();
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
        maxIter = eZmaxIter.getValue();
        //psfSplitted = eZpsfSplitted.getValue();
        double tmp;

        Sequence seqImgTmp = sequenceImg.getValue();
        Sequence seqPsfTmp = sequencePsf.getValue();
        Sequence seqResTmp = lastResult.getValue();

        if (seqImgTmp != null  && seqPsfTmp != null) {
            tmp = seqImgTmp.getSizeX(); //Just to compute the size of the coefficient, we take the width of the image
        } else {
            //Test if we have the image and the psf ...
            if(seqImgTmp == null || seqPsfTmp == null){
                //If there is a missing parameter we notify the user with the missing parameter as information
                String message = "You have forgotten to give ";
                String messageEnd = "";
                if (seqImgTmp == null) {
                    messageEnd = messageEnd.concat("the image ");
                }
                if(seqPsfTmp == null) {
                    if (seqImgTmp == null) {
                        messageEnd = messageEnd.concat("and ");
                    }
                    messageEnd = messageEnd.concat("a PSF");
                }
                message(message+messageEnd);
            }
            return;
        }
        coef = (tmp + eZcoef.getValue())/tmp;
        if (isHeadLess()) {
            reuse = false;
        }else {
            reuse = !(seqResTmp == null);
        }

        //Testing epsilon and grtol
        if (mu < 0) {
            message("Regularization level MU must be strictly positive");
        }
        if (epsilon <= 0) {
            message("Threshold level EPSILON must be strictly positive");
        }
        if (grtol < 0 || grtol >= 1) {
            message("grtol canno't be lower than 0 or greater than 1");
        }
        if (coef < 0 || coef > 3) {
            message("The Padding can not be lower than 0 or have a value greater than 3");
        }
        if (maxIter < -1)  {
            maxIter = -1;
        }
        try {
            //And if the sizes are matching
            img = seqImgTmp.getFirstNonNullImage();
            psf = seqPsfTmp.getFirstNonNullImage();
            if (reuse) {
                result = seqResTmp.getFirstNonNullImage();
            }
            //if the user they the psf is splitted and the psf and image are not of the same size
            //if (psfSplitted && (img.getWidth() != psf.getWidth() || img.getHeight() != psf.getHeight())) {
            //    message("The image and the psf should be of same size");
            //}
            //if the user make a mistake between psf and image
            if (psf.getWidth() > img.getWidth() || psf.getHeight() > img.getHeight()) {
                message("The psf can not be larger than the image");
            }
            //if the user give a bad previous result
            if (reuse && (result.getWidth() != img.getWidth() || result.getHeight() != img.getHeight() || seqResTmp.getSizeZ() != seqImgTmp.getSizeZ())) {
                message("The previous result must be the same size as the image");
            }
            //if the user does not give data with same dimensions : no 3d and 2d at same time.
            if (seqImgTmp.getSizeZ() == 1 && seqPsfTmp.getSizeZ() > 1 ||
                    seqImgTmp.getSizeZ() > 1 && seqPsfTmp.getSizeZ() == 1) {
                message("The psf and the image should have the same number of dimensions");
            }
            //if the user give data in 4D
            if (seqImgTmp.getSizeT() > 1 || seqPsfTmp.getSizeT() > 1) {
                message("Sorry we do not support 4D data for now");
            }

            //Everything seems good we are ready to launch
            if (goodInput) {
                if (reuse && tvDec != null) {
                    //If we restart, we reuse the same data and PSF
                    tvDec.setRegularizationWeight(mu);
                    tvDec.setRegularizationThreshold(epsilon);
                    tvDec.setRelativeTolerance(grtol);
                    tvDec.setMaximumIterations(maxIter);
                    tvDec.setOutputShape(shape);
                    tvDec.setPositivity(eZpositivity.getValue());
                    // We verify that the bounds are respected for previous input
                    lowerBound = tvDec.getLowerBound();
                    DoubleArray psfArray =  (DoubleArray) IcyBufferedImageUtils.imageToArray(seqPsfTmp, 0);
                    DoubleArray myArray = (DoubleArray) IcyBufferedImageUtils.imageToArray(seqResTmp, 0);
                    //DoubleArray myArray = (DoubleArray)tvDec.getResult();
                    myArray.map(new DoubleFunction() {
                        @Override
                        public double apply(double arg) {
                            if (arg >= lowerBound) {
                                return arg;
                            } else {
                                return lowerBound;
                            }
                        }
                    });
                    tvDec.setPsf(psfArray);
                    tvDec.setResult(myArray);
                    token.start();  //By default wait for the end of the job
                    computeNew = true;
                } else {
                    //Launching computation
                    tvDec = new TotalVariationJobForIcy(token);
                    tvDec.setRegularizationWeight(mu);
                    tvDec.setRegularizationThreshold(epsilon);
                    tvDec.setRelativeTolerance(grtol);
                    tvDec.setAbsoluteTolerance(gatol);
                    tvDec.setMaximumIterations(maxIter);
                    tvDec.setPositivity(eZpositivity.getValue());
                    tvDec.setViewer(new tvViewer());
                    thread.setJob(tvDec);
                    // Read the image and the PSF.
                    width = img.getWidth();
                    height = img.getHeight();
                    sizeZ = seqImgTmp.getSizeZ();
                    DoubleArray imgArray, psfArray, weight;

                    if (seqImgTmp.getSizeZ() == 1) { //2D
                        shape = Shape.make(FFTUtils.bestDimension((int)(width*coef)), FFTUtils.bestDimension((int)(height*coef)));
                    } else { //3D
                        shape = Shape.make(FFTUtils.bestDimension((int)(width*coef)),
                                FFTUtils.bestDimension((int)(height*coef)),
                                FFTUtils.bestDimension((int)(sizeZ*coef)));
                    }
                    imgArray =  (DoubleArray) IcyBufferedImageUtils.imageToArray(seqImgTmp, 0);
                    psfArray =  (DoubleArray) IcyBufferedImageUtils.imageToArray(seqPsfTmp, 0);

                    weight = createWeight(imgArray);
                    //BEWARE here we change the value to match the new padded image size
                    //addImage(weight.flatten(), "weights", width, height, sizeZ); //Uncomment this to see weights
                    width = FFTUtils.bestDimension((int)(width*coef));
                    height = FFTUtils.bestDimension((int)(height*coef));
                    sizeZ = FFTUtils.bestDimension((int)(sizeZ*coef));

                    tvDec.setWeight(weight);
                    tvDec.setData(imgArray);
                    tvDec.setPsf(psfArray);
                    tvDec.setOutputShape(shape);
                    // We verify that the bounds are respected for previous input
                    lowerBound = tvDec.getLowerBound();
                    DoubleArray myArray = tvDec.getData();
                    myArray.map(new DoubleFunction() {

                        @Override
                        public double apply(double arg) {
                            if (arg >= lowerBound) {
                                return arg;
                            } else {
                                return lowerBound;
                            }
                        }
                    });
                    token.start();  //By default wait for the end of the job
                    computeNew = true;
                }
            }
        } catch (IllegalArgumentException e) {
            new AnnounceFrame("Oops, Error: "+ e.getMessage());
        }
    }

    /****************************************************/
    /**                  UTILS FUNCTIONS               **/
    /****************************************************/
    //Small utils function that will get the sequence and convert it to ShapedArray
    private ShapedArray weightMapToArray(EzVarSequence seq){
        Sequence in = seq.getValue();
        if (in != null) {
            return IcyBufferedImageUtils.imageToArray(in, 0);
        } else {
            throw new IllegalArgumentException("The input requested was not found");
        }

    }

    /**
     * The goal is to create an weight array, but it will be created depending
     * the user input so we will have to test each cases:
     *  	-None
     *  	-A given map
     *  	-A variance map
     *  	-A computed variance
     * Then we apply the dead pixel map
     * 
     * @param data
     * @return
     */
    private DoubleArray createWeight(ShapedArray data){
        WeightGenerator weightGen = new WeightGenerator();
        String newValue = options.getValue();
        ShapedArray array;
        ShapedArray deadPixMap = null;
        if (showPixMap.getValue() && deadPixel.getValue() != null && deadPixel.getValue().getFirstNonNullImage() != null) {
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
            weightGen.setComputedVariance(data, alpha.getValue(), beta.getValue());
        } else {
            throw new IllegalArgumentException("Incorrect argument for weightmap");
        }
        weightGen.setPixelMap(deadPixMap);
        return weightGen.getWeightMap(data.getShape()).toDouble();
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

    //Copy the metadata from the input image to the output image
    //In the future if we want to change metadata should be here
    private void updateMetaData(Sequence seq) {
        Sequence imageIn = sequenceImg.getValue();
        OMEXMLMetadataImpl newMetdat = OMEUtil.createOMEMetadata(imageIn.getMetadata());
        //newMetdat.setImageDescription("MyDescription", 0);
        seq.setMetaData(newMetdat);
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
            updateMetaData(sequence);
            computeNew = false;
        }
        try {
            sequence.beginUpdate();
            double[] in = tvDec.getResult().toDouble().flatten();
            for (int j = 0; j < sizeZ; j++) {
                double[] temp = new double[width*height];
                for (int i = 0; i < width*height; i++) {
                    temp[i] = in[i+j*width*height];
                }
                sequence.setImage(0,j, new IcyBufferedImage(width, height, temp));
            }
        } finally {
            sequence.endUpdate();
        }
        sequence.setName("TV mu="+mu+" Epsilon="+epsilon+" Iteration="+tvDec.getIterations()+" grToll= "+tvDec.getRelativeTolerance());
    }

    //When the plugin is closed we try to close/stop everything
    @Override
    public void clean() {
        if (sequence != null) {
            sequence.close();
        }
        if (token != null) {
            token.stop();
            token.exit();
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
        if (token != null) {
            token.stop();
        }
    }

    //The input variable for the protocol
    @Override
    public void declareInput(VarList inputMap) {
        inputMap.add("psf", sequencePsf.getVariable());
        inputMap.add("image", sequenceImg.getVariable());
        inputMap.add("mu", eZmu.getVariable());
        inputMap.add("epsilon", eZepsilon.getVariable());
        inputMap.add("maxIter", eZmaxIter.getVariable());
        inputMap.add("coef", eZcoef.getVariable());
    }

    //The output variable for the protocol
    @Override
    public void declareOutput(VarList outputMap) {
        outputMap.add("output", output.getVariable());
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