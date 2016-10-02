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
import mitiv.array.ArrayUtils;
import mitiv.array.Double1D;
import mitiv.array.DoubleArray;
import mitiv.array.ShapedArray;
import mitiv.base.Shape;
import mitiv.base.mapping.DoubleFunction;
import mitiv.invpb.ReconstructionJob;
import mitiv.invpb.ReconstructionViewer;
import mitiv.utils.FFTUtils;
import mitiv.utils.WeightFactory;
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
import plugins.mitiv.blinddeconv.ToolTipText;
import plugins.mitiv.io.IcyBufferedImageUtils;
import plugins.mitiv.reconstruction.TotalVariationJobForIcy;
/**
 * MitivDeconvolution implements regilarized multi-dimensional deconvolution.
 */
public class MitivDeconvolution extends EzPlug implements Block, EzStoppable, SequenceListener, EzVarListener<String> {

    /****************************************************/
    /**                 VIEWER UPDATE                  **/
    /****************************************************/
    public class tvViewer implements ReconstructionViewer {
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
    private boolean validParameters = true;
    //private boolean psfSplitted = false;
    private boolean computeNew = true;
    private boolean reuse = false;

    private int width = -1;
    private int height = -1;
    private int sizeZ = -1;
    private int widthPad = -1;
    private int heightPad = -1;
    private int sizeZPad = -1;
    private double padXY = 1.0;
    private double padZ = 1.0;
    private Shape shapePad;

    private ReconstructionThreadToken token;
    ReconstructionThread thread;

    private EzVarSequence ezPsfSeq; // Point Spread Function
    private EzVarSequence ezDatSeq; // Input data
    private EzVarSequence ezResSeq; // Last result
    private EzVarSequence ezOutSeq; // Output
    
    //private EzVarBoolean eZpsfSplitted = new EzVarBoolean("Is the psf splitted?", psfSplitted);
    private EzVarDouble ezMu;
    private EzVarDouble ezEpsilon;
    private EzVarInteger ezPadXY, ezPadZ;
    private EzVarInteger ezMaxIter;
    private EzVarBoolean ezPositivity;

    private String weightOption1;
    private String weightOption2;
    private String weightOption3;
    private String weightOption4;

    private EzVarText ezWgtChoice;
    private EzVarSequence ezWgtSeq;  // statistical weights of the data
    private EzVarSequence ezVarSeq;  // variance of the data 
    private EzVarSequence ezBadSeq;  // bad data mask
    private EzVarBoolean ezShowBad;
    private EzVarDouble ezGain;      // detector gain (e-/ADU)
    private EzVarDouble ezNoise;     // standard deviation of detector noise (e-/voxel)
    private EzGroup ezRegularizationGroup;
    private EzGroup ezWeightingGroup;
    private EzGroup ezConvergenceGroup;

    //IcyBufferedImage img;
    //IcyBufferedImage psf;
    IcyBufferedImage result;

    /****************************************************/
    /**                 MESSAGE ICY                    **/
    /****************************************************/
    /**
     * Display an error message and invalidate the processing.
     * @param reason
     */
    private void error(String reason){
        new AnnounceFrame(reason);
        validParameters = false;
    }

    /****************************************************/
    /**                 INITIALIZE ICY                 **/
    /****************************************************/
    @Override
    protected void initialize() {

        weightOption1 = new String("Uniform weights");
        weightOption2 = new String("Given weights");
        weightOption3 = new String("Given variance");
        weightOption4 = new String("Computed weights");

        /* Create all widget instances. */
        ezPsfSeq     = new EzVarSequence("Input PSF");
        ezDatSeq     = new EzVarSequence("Input data");
        ezResSeq     = new EzVarSequence("Previous result");
        ezOutSeq     = new EzVarSequence("Output"); //In headLess mode only
        ezWgtSeq     = new EzVarSequence("Weights");
        ezVarSeq     = new EzVarSequence("Variance");
        ezBadSeq     = new EzVarSequence("Bad data mask");
        ezMu         = new EzVarDouble("Regularization level", 0, Double.MAX_VALUE, 0.1);
        ezEpsilon    = new EzVarDouble("L1-L2 threshold", 0, Double.MAX_VALUE, 1);
        ezPadXY      = new EzVarInteger("XY padding", 0, 10000, 1);
        ezPadZ       = new EzVarInteger("Z padding", 0, 10000, 1);
        ezMaxIter    = new EzVarInteger("Max. num. iterations", -1, Integer.MAX_VALUE, 1);
        ezPositivity = new EzVarBoolean("Enforce nonnegativity", false);
        ezWgtChoice  = new EzVarText("Options", new String[] { weightOption1, weightOption2, weightOption3, weightOption4}, 0, false);
        ezShowBad    = new EzVarBoolean("Show bad data?", false);
        ezGain       = new EzVarDouble("Gain (e-/ADU)");
        ezNoise      = new EzVarDouble("Readout noise (e-/pixel)");
        ezRegularizationGroup = new EzGroup("Regularization", ezMu, ezEpsilon);
        ezWeightingGroup      = new EzGroup("Weighting", ezWgtChoice, ezWgtSeq, ezVarSeq, ezGain, ezNoise, ezShowBad, ezBadSeq);
        ezConvergenceGroup    = new EzGroup("Convergence criterion", ezMaxIter);

        ezResSeq.setNoSequenceSelection();

        /* Set initial values. */
        ezMu.setValue(mu);
        ezEpsilon.setValue(epsilon);
        ezMaxIter.setValue(maxIter);
        ezPadXY.setValue(0);
        ezPadZ.setValue(0);
        ezWgtChoice.addVarChangeListener(this);

        /* Set visibility of weight parameters. */
        ezGain.setVisible(false);
        ezNoise.setVisible(false);
        ezBadSeq.setVisible(false);
        ezShowBad.addVisibilityTriggerTo(ezBadSeq, true);

        /* Set tool-tips. */
        ezPsfSeq.setToolTipText(ToolTipText.sequencePSF);
        ezDatSeq.setToolTipText(ToolTipText.sequenceImage);
        ezResSeq.setToolTipText(ToolTipText.booleanRestart);
        //eZpsfSplitted.setToolTipText(ToolTipText.booleanPSFSplitted);
        ezMu.setToolTipText(ToolTipText.doubleMu);
        ezEpsilon.setToolTipText(ToolTipText.doubleEpsilon);
        ezMaxIter.setToolTipText(ToolTipText.doubleMaxIter);
        ezPadXY.setToolTipText(ToolTipText.doublePadding);
        ezPadZ.setToolTipText(ToolTipText.doublePadding);
        ezGain.setToolTipText(ToolTipText.doubleGain);
        ezNoise.setToolTipText(ToolTipText.doubleNoise);
        ezPositivity.setToolTipText(ToolTipText.booleanPositivity);
        ezShowBad.setToolTipText(ToolTipText.sequencePixel);
        ezWgtSeq.setToolTipText(ToolTipText.sequenceWeigth);
        ezVarSeq.setToolTipText(ToolTipText.sequenceVariance);

        /* Add all components. */
        addEzComponent(ezDatSeq);
        addEzComponent(ezPsfSeq);
        addEzComponent(ezResSeq);
        addEzComponent(ezRegularizationGroup);
        addEzComponent(ezConvergenceGroup);
        addEzComponent(ezWeightingGroup);
        addEzComponent(ezPadXY);
        addEzComponent(ezPadZ);
        addEzComponent(ezPositivity);

        token = new ReconstructionThreadToken(new double[]{mu,epsilon,gatol,grtol});
        thread = new ReconstructionThread(token);
        thread.start();
    }

    /****************************************************/
    /**                  RUN PLUGIN                    **/
    /****************************************************/
    @Override
    protected void execute() {
        /* Extract settings. */
        validParameters = true;
        mu = ezMu.getValue();
        epsilon = ezEpsilon.getValue();
        maxIter = ezMaxIter.getValue();
        double tmp;

        Sequence datSeq = ezDatSeq.getValue();
        Sequence psfSeq = ezPsfSeq.getValue();
        Sequence resSeq = ezResSeq.getValue();

        if (datSeq == null) {
        	error("A data sequence to be deblurred must be selected");
        	return;
        }
        if (psfSeq == null) {
        	error("A point spread function (PSF) must be selected");
        	return;
        }
        tmp = datSeq.getSizeX(); //Just to compute the size of the coefficient, we take the width of the image
        padXY= (tmp + ezPadXY.getValue())/tmp;
        padZ = (datSeq.getSizeZ() + ezPadZ.getValue())/datSeq.getSizeZ();
        if (isHeadLess()) {
            reuse = false;
        } else {
            reuse = (resSeq != null);
        }

        //Testing epsilon and grtol
        if (mu < 0) {
            error("Regularization level MU must be strictly positive");
        }
        if (epsilon <= 0) {
            error("Threshold level EPSILON must be strictly positive");
        }
        if (grtol < 0 || grtol >= 1) {
            error("grtol canno't be lower than 0 or greater than 1");
        }
        if (padXY < 0|| padZ < 0) {
            error("The Padding can not be lower than 0");
        }
        if (maxIter < 0)  {
            maxIter = -1;
        }
        try {
            //And if the sizes are matching
            if (reuse) {
                result = resSeq.getFirstNonNullImage();
            }
            //if the user they the psf is splitted and the psf and image are not of the same size
            //if (psfSplitted && (img.getWidth() != psf.getWidth() || img.getHeight() != psf.getHeight())) {
            //    message("The image and the psf should be of same size");
            //}
            //if the user make a mistake between psf and image
            if (psfSeq.getWidth() > datSeq.getWidth() || psfSeq.getHeight() > datSeq.getHeight()) {
                error("The psf can not be larger than the image");
            }
            //if the user give a bad previous result
            if(reuse) {
                boolean sameAsOrigin = resSeq.getSizeX() == datSeq.getSizeX() 
                        && resSeq.getSizeY() == datSeq.getSizeY() 
                        && resSeq.getSizeZ() == datSeq.getSizeZ();
                boolean sameAsPrevious = resSeq.getSizeX() == shapePad.dimension(0) 
                        && resSeq.getSizeY() == shapePad.dimension(1) 
                        && (resSeq.getSizeZ() == 1 ? true : resSeq.getSizeZ() == shapePad.dimension(2)); // If we are a 2d image, we do nothing
                if (!(sameAsOrigin || sameAsPrevious)) {
                    error("The previous result does not have the same dimensions as the input image");
                }
            }
            //if the user does not give data with same dimensions : no 3d and 2d at same time.
            if (datSeq.getSizeZ() == 1 && psfSeq.getSizeZ() > 1 ||
                    datSeq.getSizeZ() > 1 && psfSeq.getSizeZ() == 1) {
                error("The psf and the image should have the same number of dimensions in Z");
            }
            //if the user give data in 4D
            if (datSeq.getSizeT() > 1 || psfSeq.getSizeT() > 1) {
                error("Sorry we do not support 4D data for now");
            }

            //Everything seems good we are ready to launch
            if (validParameters) {
                // Settings all sizes
                width  = Math.max(datSeq.getWidth(), psfSeq.getWidth());
                height = Math.max(datSeq.getHeight(), psfSeq.getHeight());
                sizeZ  = Math.max(datSeq.getSizeZ(), psfSeq.getSizeZ());
                widthPad  = FFTUtils.bestDimension((int)(width*padXY));
                heightPad = FFTUtils.bestDimension((int)(height*padXY));
                sizeZPad  = FFTUtils.bestDimension((int)(sizeZ*padZ));

                if (datSeq.getSizeZ() == 1) { //2D
                    // shape = Shape.make(width, height);
                    shapePad = Shape.make(widthPad, heightPad);
                } else { //3D
                    // shape = Shape.make(width, height, sizeZ);
                    shapePad = Shape.make(widthPad, heightPad, sizeZPad);
                }
                // FIXME What happen when the image is not square, should we pad it ?
                if (reuse && tvDec != null) {
                    //If we restart, we reuse the same data and PSF
                    tvDec.setRegularizationWeight(mu);
                    tvDec.setRegularizationThreshold(epsilon);
                    tvDec.setRelativeTolerance(grtol);
                    tvDec.setMaximumIterations(maxIter);
                    tvDec.setOutputShape(shapePad);
                    tvDec.setPositivity(ezPositivity.getValue());
                    // We verify that the bounds are respected for previous input
                    lowerBound = tvDec.getLowerBound();
                    DoubleArray psfArray =  (DoubleArray) IcyBufferedImageUtils.imageToArray(psfSeq, 0);
                    DoubleArray resArray = (DoubleArray) IcyBufferedImageUtils.imageToArray(resSeq, 0);
                    //DoubleArray myArray = (DoubleArray)tvDec.getResult();
                    resArray.map(new DoubleFunction() {
                        @Override
                        public double apply(double arg) {
                            if (arg >= lowerBound) {
                                return arg;
                            } else {
                                return lowerBound;
                            }
                        }
                    });
                    resArray = (DoubleArray) ArrayUtils.pad(resArray, shapePad);
                    tvDec.setPsf(psfArray);
                    tvDec.setResult(resArray);
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
                    tvDec.setPositivity(ezPositivity.getValue());
                    tvDec.setViewer(new tvViewer());
                    thread.setJob(tvDec);
                    // Read the image and the PSF.

                    DoubleArray imgArray, psfArray, weight;

                    imgArray =  (DoubleArray) IcyBufferedImageUtils.imageToArray(datSeq, 0);
                    psfArray =  (DoubleArray) IcyBufferedImageUtils.imageToArray(psfSeq, 0);
                    weight = createWeight(imgArray);

                    imgArray = (DoubleArray) ArrayUtils.pad(imgArray, shapePad);
                    psfArray = (DoubleArray) ArrayUtils.pad(psfArray, shapePad);
                    weight   = (DoubleArray) ArrayUtils.pad(weight  , shapePad);

                    //BEWARE here we change the value to match the new padded image size
                    //addImage(weight.flatten(), "weights", widthPad, heightPad, sizeZPad); //Uncomment this to see weights

                    tvDec.setWeight(weight);
                    tvDec.setData(imgArray);
                    tvDec.setPsf(psfArray);
                    tvDec.setOutputShape(shapePad);
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
            e.printStackTrace();
            new AnnounceFrame("Oops, Error: "+ e.getMessage());
        }
    }

    /****************************************************/
    /**                  UTILS FUNCTIONS               **/
    /****************************************************/
    //Small utils function that will get the sequence and convert it to ShapedArray
    private static ShapedArray weightMapToArray(EzVarSequence seq){
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
    private DoubleArray createWeight(ShapedArray datArray) {
    	return null;
//        String option = options.getValue();
//    	ShapedArray wgtArray = null;
//        Sequence seq;
//        boolean wgtCopy = true, datCopy = true;
//
//        ShapedArray array;
//        ShapedArray deadPixMap = null;
//        if (option == weightOption1) {
//        	//No options
//            //First case, we create an array of 1, that correspond to the image
//            double[] weight = new double[width*height*sizeZ];
//            for (int i = 0; i < weight.length; i++) {
//                weight[i] = 1;
//            }
//            weightGen.setWeightMap(Double1D.wrap(weight, weight.length));
//        } else if (option == weightOption2) {//Weight Map
//            //If the user give a weight map we convert it and give it to weightGen
//            array = weightMapToArray(weightMap);
//            weightGen.setWeightMap(array);
//        } else if(option == weightOption3) {//Variance map
//            //If the user give a varianceMap map we give it to weightGen
//            array = weightMapToArray(varianceMap);
//            weightGen.setWeightMap(array);//Gain + readOut
//        } else if(option == weightOption4) {
//            weightGen.setComputedVariance(data, gain.getValue(), noise.getValue());
//        } else {
//            throw new IllegalArgumentException("Incorrect argument for weightmap");
//        }
//        weightGen.setPixelMap(deadPixMap);
//        return weightGen.getWeightMap(data.getShape()).toDouble();
//        
//
//        // We check the values given
//        if (options.getValue() == weightOption1) {
//        	// Nothing specified.  The weights are an array of ones with same size as the data.
//        	wgtArray = WeightFactory.defaultWeights(datArray);
//        	wgtCopy = false; // no needs to copy weights
//        } else if (options.getValue() == weightOption2) {
//        	// A map of weights is provided.
//            if ((seq = weights.getValue()) != null) {      
//                wgtArray = IcyBufferedImageUtils.imageToArray(seq.getAllImage());
//            }
//        } else if (options.getValue() == weightOption3) {
//        	// A variance map is provided. FIXME: check shape and values.
//            if ((seq = weights.getValue()) != null) {
//            	ShapedArray varArray = IcyBufferedImageUtils.imageToArray(weightMap.getAllImage());
//                wgtArray = WeightFactory.computeWeightsFromVariance(varArray);
//            	wgtCopy = false; // no needs to copy weights
//            }
//        } else if (options.getValue() == weightOption4) {
//        	// Weights are computed given the gain and the readout noise of the detector.
//        	double gamma = gain.getValue();
//        	double sigma = noise.getValue();
//        	double alpha = 1/gamma;
//        	double beta = (sigma/gamma)*(sigma/gamma);
//            wgtArray = WeightFactory.computeWeightsFromData(datArray, alpha, beta);
//        	wgtCopy = false; // no needs to copy weights
//        }
//        
//        // Make sure weights and data are private copies because we may have to modify their contents.
//        if (wgtCopy) {
//        	wgtArray = flatCopy(wgtArray);
//        }
//        if (datCopy) {
//        	datArray = flatCopy(datArray);
//        }
//        
//        if (showPixMap.getValue() && deadPixel.getValue() != null && deadPixel.getValue().getFirstNonNullImage() != null) {
//            deadPixMap = weightMapToArray(deadPixel);
//        }
//        if (deadPixGiven.getValue() && (seq = deadPixel.getValue()) != null) {
//        	// Account for bad data.
//        	ShapedArray badArr = IcyBufferedImageUtils.imageToArray(seq.getAllImage());
//        	WeightFactory.removeBads(wgtArray, badArr);
//        }
//
//    	// Check everything.
//    	WeightFactory.fixWeightsAndData(wgtArray, datArray);
//    	return wgtArray;

        
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
        Sequence imageIn = ezDatSeq.getValue();
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
                ezOutSeq.setValue(sequence);
            }else{
                addSequence(sequence);
            }
            updateMetaData(sequence);
            computeNew = false;
        }
        try {
            sequence.beginUpdate();
            double[] in = tvDec.getResult().toDouble().flatten();
            for (int j = 0; j < sizeZPad; j++) {
                double[] temp = new double[widthPad*heightPad];
                for (int i = 0; i < widthPad*heightPad; i++) {
                    temp[i] = in[i+j*widthPad*heightPad];
                }
                sequence.setImage(0,j, new IcyBufferedImage(widthPad, heightPad, temp));
            }
        } finally {
            sequence.endUpdate();
        }
        sequence.setName("TV mu="+mu+" Epsilon="+epsilon+" Iteration="+tvDec.getIterations());
    }

    //When the plugin is closed we try to close/stop everything
    @Override
    public void clean() {
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
            ezWgtSeq.setVisible(false);
            ezVarSeq.setVisible(false);
            ezGain.setVisible(false);
            ezNoise.setVisible(false);
        } else if(newValue == weightOption2) {
            ezWgtSeq.setVisible(true);
            ezVarSeq.setVisible(false);
            ezGain.setVisible(false);
            ezNoise.setVisible(false);
        } else if(newValue == weightOption3) {
            ezWgtSeq.setVisible(false);
            ezVarSeq.setVisible(true);
            ezGain.setVisible(false);
            ezNoise.setVisible(false);
        } else if(newValue == weightOption4) {
            ezWgtSeq.setVisible(false);
            ezVarSeq.setVisible(false);
            ezGain.setVisible(true);
            ezNoise.setVisible(true);
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
        initialize();
        inputMap.add("image", ezDatSeq.getVariable());
        inputMap.add("psf", ezPsfSeq.getVariable());
        inputMap.add("mu", ezMu.getVariable());
        inputMap.add("epsilon", ezEpsilon.getVariable());
        inputMap.add("maxIter", ezMaxIter.getVariable());
        inputMap.add("coefXY", ezPadXY.getVariable());
        inputMap.add("coefZ", ezPadZ.getVariable());
    }

    //The output variable for the protocol
    @Override
    public void declareOutput(VarList outputMap) {
        outputMap.add("output", ezOutSeq.getVariable());
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