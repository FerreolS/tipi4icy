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

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.SwingUtilities;

import icy.gui.frame.progress.AnnounceFrame;
import icy.image.IcyBufferedImage;
import icy.sequence.Sequence;
import mitiv.array.Array3D;
import mitiv.array.ArrayFactory;
import mitiv.array.ArrayUtils;
import mitiv.array.DoubleArray;
import mitiv.array.FloatArray;
import mitiv.array.ShapedArray;
import mitiv.base.Shape;
import mitiv.base.Traits;
import mitiv.invpb.EdgePreservingDeconvolution;
import mitiv.linalg.shaped.ShapedVector;
import mitiv.optim.OptimTask;
import mitiv.utils.FFTUtils;
import mitiv.utils.WeightFactory;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzButton;
import plugins.adufour.ezplug.EzGroup;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;
import plugins.adufour.ezplug.EzVar;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarChannel;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarDoubleArrayNative;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarListener;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.ezplug.EzVarText;
import plugins.mitiv.io.IcyBufferedImageUtils;
/**
 * MitivDeconvolution implements regularized multi-dimensional deconvolution.
 */
public class MitivDeconvolution extends EzPlug implements Block, EzStoppable {


    /***************************************************/
    /**         Plugin interface variables            **/
    /***************************************************/
    private EzVarSequence image;
    private EzVarChannel channel;
    private EzVarSequence psf; // Point Spread Function

    private EzVarInteger    paddingSizeXY, paddingSizeZ;
    private EzVarText imageSize, outputSize;
    /** weighting tab: **/
    private EzVarText weightsMethod;
    private final String[] weightOptions = new String[]{"None","Inverse covariance map","Variance map","Computed variance"};
    private EzVarDouble  gain, noise;
    private EzVarSequence weights, deadPixel;
    private EzButton showWeight;


    /** deconvolution tab: **/
    private EzVarDouble logmu, mu, epsilon;
    private EzVarSequence  restart;
    private EzVarBoolean  positivity;
    private EzVarBoolean  singlePrecision;
    private EzVarBoolean  showIteration;
    private EzButton startDec, stopDec, showFull;
    private EzVarInteger    nbIterDeconv;

    /** headless mode: **/
    private EzVarSequence   outputHeadlessImage;



    /****************************************************/
    /**                 VARIABLES                      **/
    /****************************************************/

    static boolean debug =false;

    private int sizeX=512, sizeY=512, sizeZ=128; // Input sequence size
    protected int psfSizeX=1,psfSizeY=1,psfSizeZ=1;
    private Shape imageShape, psfShape;
    private  int Nxy=512, Nz=128;            // Output (padded sequence size)
    private Shape outputShape;
    private static double[][] scaleDef =new double[][] {{1.0},{1.0 ,1.0},{1.0 ,1.0, 1.0},{1.0 ,1.0, 1.0,1.0}};

    Sequence cursequence=null;
    EdgePreservingDeconvolution solver =  new EdgePreservingDeconvolution();
    private boolean run;
    private EzVarChannel channelpsf, channelRestart;
    private EzGroup ezPaddingGroup;
    private EzGroup ezWeightingGroup;
    private EzGroup ezDeconvolutionGroup;
    private EzVarDoubleArrayNative scale;
    private EzGroup ezDeconvolutionGroup2;


    protected Sequence lastSequence;


    /*********************************/
    /**      Initialization         **/
    /*********************************/
    @Override
    protected void initialize() {
        if (!isHeadLess()) {
            getUI().setParametersIOVisible(false);
            getUI().setActionPanelVisible(false);
            outputHeadlessImage = new EzVarSequence("Output Image");
        }


        image = new EzVarSequence("Image:");
        channel = new EzVarChannel("Image channel:", image.getVariable(), false);
        psf = new EzVarSequence("PSF:");
        channelpsf = new EzVarChannel("PSF channel :", psf.getVariable(), false);
        restart = new EzVarSequence("Starting point:");
        restart.setNoSequenceSelection();
        channelRestart = new EzVarChannel("Initialization channel :", restart.getVariable(), false);

        imageSize = new EzVarText("Image size:");
        outputSize = new EzVarText("Output size:");
        paddingSizeXY = new EzVarInteger("padding xy:",0, Integer.MAX_VALUE,1);
        paddingSizeZ = new EzVarInteger("padding z :",0, Integer.MAX_VALUE,1);

        image.setNoSequenceSelection();
        psf.setNoSequenceSelection();
        ezPaddingGroup = new EzGroup("Padding",paddingSizeXY,paddingSizeZ);
        ezPaddingGroup.setFoldedState(true);

        EzVarListener<Integer> zeroPadActionListener = new EzVarListener<Integer>() {
            @Override
            public void variableChanged(EzVar<Integer> source, Integer newValue) {
                updatePaddedSize();
                updateImageSize();
                updateOutputSize();
            }
        };
        paddingSizeXY.addVarChangeListener(zeroPadActionListener);
        paddingSizeZ.addVarChangeListener(zeroPadActionListener);




        restart.addVarChangeListener(new EzVarListener<Sequence>() {
            @Override
            public void variableChanged(EzVar<Sequence> source,
                    Sequence newValue) {
                newValue = restart.getValue();
                if (newValue != null || (newValue != null && newValue.isEmpty())) {
                    System.out.println("restart changed:"+newValue.getName());
                }
            }
        });
        image.addVarChangeListener(new EzVarListener<Sequence>() {
            @Override
            public void variableChanged(EzVar<Sequence> source,
                    Sequence newValue) {
                if (debug) {
                    System.out.println("Seq ch..."+image.getValue());
                }

                Sequence seq = image.getValue();
                if (seq != null || (seq != null && seq.isEmpty())) {
                    sizeX =  newValue.getSizeX();
                    sizeY = newValue.getSizeY();
                    sizeZ = newValue.getSizeZ();
                    updatePSFSize();
                    updateImageSize();

                    imageShape = new Shape(sizeX, sizeY, sizeZ);


                    if (debug) {
                        System.out.println("Seq changed:" + sizeX + "  "+ Nxy);
                    }
                    // setting restart value to the current sequence
                    restart.setValue(newValue);
                    channelRestart.setValue(channel.getValue());
                }
            }

        });

        channel.addVarChangeListener(new EzVarListener<Integer>() {
            @Override
            public void variableChanged(EzVar<Integer> source, Integer newValue) {
                channelRestart.setValue(newValue);
            }
        });

        psf.addVarChangeListener(new EzVarListener<Sequence>() {
            @Override
            public void variableChanged(EzVar<Sequence> source,
                    Sequence newValue) {
                if (debug) {
                    System.out.println("PSF changed"+psf.getValue());
                }
                if (newValue != null || (newValue != null && newValue.isEmpty())) {
                    psfSizeX = Math.max(1,newValue.getSizeX());
                    psfSizeY =  Math.max(1,newValue.getSizeY());
                    psfSizeZ =  Math.max(1,newValue.getSizeZ());
                    psfShape = new Shape(psfSizeX, psfSizeY, psfSizeZ);
                    updatePSFSize();
                    updateImageSize();

                    if (debug) {
                        System.out.println("PSF changed:" + psfSizeX + "  "+ psfSizeY);
                    }
                }
            }

        });






        /****************************************************/
        /**                WEIGHTING GROUP                   **/
        /****************************************************/
        weightsMethod = new EzVarText(      "Weighting:", weightOptions,0, false);
        weights = new EzVarSequence(        "Map:");
        gain = new EzVarDouble(             "Gain:",1.,0.01,Double.MAX_VALUE,1);
        noise = new EzVarDouble(            "Readout Noise:",10.,0.,Double.MAX_VALUE,0.1);
        deadPixel = new EzVarSequence(      "Bad data map:");
        weights.setNoSequenceSelection();

        weightsMethod.addVarChangeListener(new EzVarListener<String>() {
            @Override
            public void variableChanged(EzVar<String> source, String newValue) {
                if (weightsMethod.getValue() == weightOptions[0]) { //None
                    weights.setVisible(false);
                    gain.setVisible(false);
                    noise.setVisible(false);
                } else if (weightsMethod.getValue() == weightOptions[1] || weightsMethod.getValue() == weightOptions[2]) {  //Personnalized map ou Variance map
                    weights.setVisible(true);
                    gain.setVisible(false);
                    noise.setVisible(false);
                    weights.setNoSequenceSelection();
                } else if (weightsMethod.getValue() == weightOptions[3]) {  //Computed variance
                    weights.setVisible(false);
                    gain.setVisible(true);
                    noise.setVisible(true);
                    weights.setNoSequenceSelection();
                } else {
                    throwError("Invalid argument passed to weight method");
                    return;
                }
            }
        });


        weights.setVisible(false);
        gain.setVisible(false);
        noise.setVisible(false);
        deadPixel.setVisible(true);
        deadPixel.setNoSequenceSelection();
        showWeight = new EzButton("Show weight map", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                DoubleArray imgArray, wgtArray;
                // Preparing parameters and testing input
                Sequence imgSeq = image.getValue();
                imgArray = (DoubleArray) IcyBufferedImageUtils.imageToArray(imgSeq, imageShape, channel.getValue());
                wgtArray = createWeights(imgArray).toDouble();
                show(wgtArray,"Weight map");
                if (debug) {
                    System.out.println("Weight compute");

                }
            }
        });
        ezWeightingGroup = new EzGroup("Weighting",weightsMethod,weights,gain,noise,deadPixel,showWeight);
        ezWeightingGroup.setFoldedState(true);



        /****************************************************/
        /**                    DECONV GROUP                  **/
        /****************************************************/
        mu = new EzVarDouble("Regularization level:",1E-5,0.,Double.MAX_VALUE,0.01);
        logmu = new EzVarDouble("Log10 of the Regularization level:",-5,-Double.MAX_VALUE,Double.MAX_VALUE,1);
        epsilon = new EzVarDouble("Threshold level:",1E-2,0.,Double.MAX_VALUE,0.01);
        nbIterDeconv = new EzVarInteger("Number of iterations: ",10,1,Integer.MAX_VALUE ,1);
        positivity = new EzVarBoolean("Enforce nonnegativity:", true);
        singlePrecision = new EzVarBoolean("Compute in single precision:", false);
        showIteration = new EzVarBoolean("Show intermediate results:", true);

        if (isHeadLess()){
            showIteration.setValue(false);
        }

        scale = new EzVarDoubleArrayNative("Aspect ratio of a voxel", scaleDef, 2,true); //FIXME SCALE
        ezDeconvolutionGroup2 = new EzGroup("More  parameters",epsilon,scale,positivity,singlePrecision);
        ezDeconvolutionGroup2.setFoldedState(true);
        ezDeconvolutionGroup = new EzGroup("Deconvolution parameters",logmu,mu,nbIterDeconv,ezDeconvolutionGroup2);

        startDec = new EzButton("Start Deconvolution", new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                Thread workerThread = new Thread() {
                    @Override
                    public void run() {
                        launch();
                    }
                };
                workerThread.start();
            }
        });
        stopDec = new EzButton("STOP Computation", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                stopExecution();
            }
        });

        showFull = new EzButton("Show the full (padded) object", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(debug){
                    System.out.println("showFull");
                }
                Sequence fSeq;
                fSeq = new Sequence("Deconvolved image");
                fSeq.copyMetaDataFrom(image.getValue(), false);
                show(solver.getSolution(),fSeq,"Deconvolved image: mu="+solver.getRegularizationLevel() );
            }
        });

        EzVarListener<Double> logmuActionListener = new EzVarListener<Double>() {
            @Override
            public void variableChanged(EzVar<Double> source, Double newValue) {
                mu.setValue(Math.pow(10, logmu.getValue()));
            }
        };
        logmu.addVarChangeListener(logmuActionListener);

        EzGroup groupVisu  = new EzGroup("Visualization", showIteration,showFull);
        groupVisu.setFoldedState(true);

        EzGroup groupStop1 = new EzGroup("Emergency STOP", stopDec);

        addEzComponent(image);
        addEzComponent(channel);
        addEzComponent(psf);
        addEzComponent(channelpsf);
        addEzComponent(restart);
        addEzComponent(channelRestart);
        addEzComponent(imageSize);
        addEzComponent(outputSize);
        addEzComponent(ezPaddingGroup);
        addEzComponent(ezWeightingGroup);
        addEzComponent(ezDeconvolutionGroup);
        addEzComponent(startDec);
        addEzComponent(groupVisu);
        addEzComponent(groupStop1);


        outputSize.setEnabled(false);
        imageSize.setEnabled(false);
        mu.setEnabled(false);


    }

    protected  void show(ShapedVector  arr) {
        show(  arr.asShapedArray(),  null,  "" ) ;
    }
    protected  void show(ShapedVector  arr,  String title ) {
        show(  arr.asShapedArray(),  null,  title ) ;
    }
    protected  void show(ShapedArray  arr,  String title ) {
        show(  arr,  null,  title ) ;
    }
    protected void show(ShapedVector  arr, Sequence sequence, String title ) {
        show(  arr.asShapedArray(),  sequence,  title ) ;
    }
    protected  void show(ShapedArray  arr) {
        show(  arr,  null,  "" ) ;
    }
    protected void show(ShapedArray  arr, Sequence sequence, String title ) {

        if (sequence == null )  {

            sequence = new Sequence();
            if (!isHeadLess()){
                addSequence(sequence);
            }
        }

        if( sequence.getFirstViewer() == null){
            if (!isHeadLess()){
                addSequence(sequence);
            }
        }
        sequence.beginUpdate();

        switch (arr.getRank()) {
            case 2:
                sequence.setImage(0,0, new IcyBufferedImage(arr.getDimension(0), arr.getDimension(1), arr.flatten()));
                break;
            case 3:
                for (int j = 0; j < arr.getDimension(2); j++) {
                    sequence.setImage(0,j, new IcyBufferedImage(arr.getDimension(0), arr.getDimension(1),((Array3D)arr).slice(j).flatten() ));
                }
                break;

            case 4:

                for (int k = 0; k < arr.getDimension(3); k++) {
                    for (int j = 0; j < arr.getDimension(2); j++) {
                        sequence.setImage(k,j, new IcyBufferedImage(arr.getDimension(0), arr.getDimension(1),((Array3D)arr).slice(k).slice(j).flatten() ));
                    }
                }
            default:
                throwError("Show can plot only 2D to 4D arrays");
                break;
        }

        sequence.endUpdate();
        sequence.setName(title);



    }

    protected void updateOutputSize() {
        String text = Nxy+"x"+Nxy+"x"+Nz;
        outputSize.setValue(text);
    }

    protected void updateImageSize() {
        String text = sizeX+"x"+sizeY+"x"+sizeZ;
        imageSize.setValue(text);
    }


    protected void updatePaddedSize() {
        int sizeXY = Math.max(sizeX, sizeY);
        Nxy = FFTUtils.bestDimension(sizeXY + paddingSizeXY.getValue());
        Nz= FFTUtils.bestDimension(sizeZ + paddingSizeZ.getValue());
        outputShape = new Shape(Nxy, Nxy, Nz);
        updateOutputSize();
        if(debug){
            System.out.println(" UpdatePaddedSize" + paddingSizeXY.getValue()  + outputShape.toString());
        }
    }

    protected void updatePSFSize() {
        paddingSizeXY.setValue( Math.max(psfSizeX, psfSizeY));
        paddingSizeZ.setValue(  psfSizeZ);
        updatePaddedSize();
        if(debug){
            System.out.println(" UpdatePaddedSize " + paddingSizeXY.getValue()  + outputShape.toString());
        }
    }

    /****************************************************/
    /**                  RUN PLUGIN                    **/
    /****************************************************/
    @Override
    protected void execute() {
        launch();
    }

    protected void launch() {
        solver = new EdgePreservingDeconvolution();

        // Preparing parameters and testing input
        Sequence imgSeq = image.getValue();
        Sequence psfSeq = psf.getValue();
        Sequence restartSeq = restart.getValue();

        psfShape = new Shape(psfSizeX, psfSizeY, psfSizeZ);


        if (imgSeq == null)
        {
            throwError("An image should be given");
            return;
        }
        if (psfSeq == null)
        {
            throwError("A psf should be given");
            return;
        }

        // Set the informations about the input
        if (sizeZ == 1) {
            throwError("Input data must be 2D or 3D");
            return;
        }
        if (paddingSizeXY.getValue() < 0.0) {
            throwError("Padding value cannot be negative");
            return;
        }
        if (paddingSizeZ.getValue() < 0.0) {
            throwError("Padding value cannot be negative");
            return;
        }

        ShapedArray imgArray, psfArray, wgtArray,restartArray;
        imgArray = IcyBufferedImageUtils.imageToArray(imgSeq, imageShape, channel.getValue());
        psfArray = IcyBufferedImageUtils.imageToArray(psfSeq, psfShape, channelpsf.getValue());

        if (restart.getValue() != null && restartSeq != null){
            restartArray = IcyBufferedImageUtils.imageToArray(restartSeq, imageShape, channelRestart.getValue());
            System.out.println("restart seq:" +restartSeq.getName());
            solver.setInitialSolution(restartArray);
        }else{
            System.out.println("restart seq is null:");
        }

        if(singlePrecision.getValue()){
            wgtArray = createWeights(imgArray).toFloat();
            solver.setForceSinglePrecision(true);
        }else{
            wgtArray = createWeights(imgArray).toDouble();
            solver.setForceSinglePrecision(false);
        }
        cursequence = new Sequence("Current Iterate");
        cursequence.copyMetaDataFrom(imgSeq, false);

        solver.setRelativeTolerance(0.0);
        solver.setUseNewCode(false);
        solver.setObjectShape(outputShape);
        solver.setPSF(psfArray);
        solver.setData(imgArray);
        solver.setWeight(wgtArray);
        solver.setEdgeThreshold(epsilon.getValue());
        solver.setRegularizationLevel(mu.getValue());
        if (scale.getValue().length !=3){
            throwError("Pixel scale must have 3 elements");
            return;
        }

        solver.setScale(scale.getValue());
        solver.setSaveBest(true);
        solver.setLowerBound(positivity.getValue() ? 0.0 : Double.NEGATIVE_INFINITY);
        solver.setUpperBound(Double.POSITIVE_INFINITY);
        solver.setMaximumIterations(nbIterDeconv.getValue());
        solver.setMaximumEvaluations(2*nbIterDeconv.getValue());
        System.out.println("Launch it:"+nbIterDeconv.getValue());
        run = true;
        OptimTask task = solver.start();

        while (run) {
            task = solver.getTask();
            System.out.println("  it "+solver.getIterations());
            if (task == OptimTask.ERROR) {
                System.err.format("Error: %s\n", solver.getReason());
                break;
            }
            if (task == OptimTask.NEW_X || task == OptimTask.FINAL_X) {
                if(showIteration.getValue()){
                    show(ArrayUtils.crop(solver.getSolution(),imageShape),cursequence,"Current mu="+solver.getRegularizationLevel() +"it:"+solver.getIterations());
                }
                if (task == OptimTask.FINAL_X) {
                    break;
                }
            }
            if (task == OptimTask.WARNING) {
                break;
            }
            solver.iterate();
        }
        show(ArrayUtils.crop(solver.getBestSolution().asShapedArray(),imageShape),cursequence,"Deconvolved image,  mu="+solver.getRegularizationLevel());

        if (isHeadLess()) {
            outputHeadlessImage.setValue(cursequence);
        }

        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                if(debug){
                    System.out.println("invoke later");
                }
                restart.setValue(cursequence);
                channelRestart.setValue(0);
            }
        });
    }

    /****************************************************/
    /**                  UTILS FUNCTIONS               **/
    /****************************************************/
    private static void throwError(String s){
        new AnnounceFrame(s);
        //throw new IllegalArgumentException(s);
    }

    /**
     * The goal is to create an array of weights, but it will be created depending
     * the user input so we will have to test each cases:
     *      - None
     *      - A given map
     *      - A variance map
     *      - A computed variance
     * Then we apply the dead pixel map
     *
     * @param datArray - The data to deconvolve.
     * @return The weights.
     */
    private ShapedArray createWeights(ShapedArray datArray) {
        ShapedArray wgtArray = null;
        Sequence seq;
        boolean wgtCopy = true, datCopy = true;

        // We check the values given
        if (weightsMethod.getValue() == weightOptions[0]) {
            // Nothing specified.  The weights are an array of ones with same size as the data.
            wgtArray = WeightFactory.defaultWeights(datArray);
            wgtCopy = false; // no needs to copy weights
        } else if (weightsMethod.getValue() == weightOptions[1]) {
            // A map of weights is provided.
            if ((seq = weights.getValue()) != null) {
                wgtArray = IcyBufferedImageUtils.imageToArray(seq.getAllImage());
            }
        } else if (weightsMethod.getValue() == weightOptions[2]) {
            // A variance map is provided. FIXME: check shape and values.
            if ((seq = weights.getValue()) != null) {
                ShapedArray varArray = IcyBufferedImageUtils.imageToArray(seq.getAllImage());
                wgtArray = WeightFactory.computeWeightsFromVariance(varArray);
                wgtCopy = false; // no needs to copy weights
            }
        } else if (weightsMethod.getValue() == weightOptions[3]) {
            // Weights are computed given the gain and the readout noise of the detector.
            double gamma = gain.getValue();
            double sigma = noise.getValue();
            double alpha = 1/gamma;
            double beta = (sigma/gamma)*(sigma/gamma);
            wgtArray = WeightFactory.computeWeightsFromData(datArray, alpha, beta);
            wgtCopy = false; // no needs to copy weights
        }

        // Make sure weights and data are private copies because we may have to modify their contents.
        if (wgtCopy) {
            wgtArray = flatCopy(wgtArray);
        }
        if (datCopy) {
            datArray = flatCopy(datArray);
        }

        if (/*deadPixGiven.getValue() && */(seq = deadPixel.getValue()) != null) {
            // Account for bad data.
            ShapedArray badArr = IcyBufferedImageUtils.imageToArray(seq.getAllImage());
            WeightFactory.removeBads(wgtArray, badArr);
        }

        // Check everything.
        WeightFactory.fixWeightsAndData(wgtArray, datArray);
        return wgtArray;
    }

    /**
     * Make a private flat copy of an array.
     *
     * @param arr - The source array.
     *
     * @return A copy of the source array.
     */
    private static ShapedArray flatCopy(ShapedArray arr)
    {
        switch (arr.getType()) {
            case Traits.FLOAT:
                return ArrayFactory.wrap(((FloatArray)arr).flatten(true), arr.getShape());
            case Traits.DOUBLE:
                return ArrayFactory.wrap(((DoubleArray)arr).flatten(true), arr.getShape());
            default:
                throw new IllegalArgumentException("Unsupported data type");
        }
    }



    //If the user call the stop button
    @Override
    public void stopExecution() {
        run = false;

    }

    //The input variable for the protocol
    @Override
    public void declareInput(VarList inputMap) {
        initialize();
        inputMap.add("image", image.getVariable());
        inputMap.add("image channel", channel.getVariable());
        inputMap.add("psf", psf.getVariable());
        inputMap.add("psf channel", channelpsf.getVariable());
        inputMap.add("starting point", restart.getVariable());
        inputMap.add("starting point channel", channelRestart.getVariable());

        inputMap.add("Padding XY", paddingSizeXY.getVariable());
        inputMap.add("Padding Z", paddingSizeZ.getVariable());

        inputMap.add("deadPixel", deadPixel.getVariable());
        inputMap.add("gain", gain.getVariable());
        inputMap.add("noise", noise.getVariable());

        inputMap.add("mu", mu.getVariable());
        inputMap.add("espilon", epsilon.getVariable());
        inputMap.add("Postivity", positivity.getVariable());
        inputMap.add("nbIteration", nbIterDeconv.getVariable());
        inputMap.add("positivity", positivity.getVariable());
        inputMap.add("single precision", singlePrecision.getVariable());
    }

    //The output variable for the protocol
    @Override
    public void declareOutput(VarList outputMap) {
        outputMap.add("outputSize", outputSize.getVariable());
        outputMap.add("output", outputHeadlessImage.getVariable());
    }

    @Override
    public void clean() {
        // TODO Auto-generated method stub

    }
}