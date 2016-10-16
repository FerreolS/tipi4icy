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

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import icy.gui.frame.progress.AnnounceFrame;
import icy.image.IcyBufferedImage;
import icy.sequence.Sequence;
import mitiv.array.Array3D;
import mitiv.array.ArrayFactory;
import mitiv.array.DoubleArray;
import mitiv.array.FloatArray;
import mitiv.array.ShapedArray;
import mitiv.base.Shape;
import mitiv.base.Traits;
import mitiv.invpb.EdgePreservingDeconvolution;
import mitiv.linalg.Vector;
import mitiv.linalg.shaped.ShapedVector;
import mitiv.optim.OptimTask;
import mitiv.utils.FFTUtils;
import mitiv.utils.WeightFactory;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzButton;
import plugins.adufour.ezplug.EzGroup;
import plugins.adufour.ezplug.EzLabel;
import plugins.adufour.ezplug.EzPanel;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;
import plugins.adufour.ezplug.EzTabs;
import plugins.adufour.ezplug.EzTabs.TabPlacement;
import plugins.adufour.ezplug.EzVar;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarChannel;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarListener;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.ezplug.EzVarText;
import plugins.mitiv.io.IcyBufferedImageUtils;
/**
 * MitivDeconvolution implements regilarized multi-dimensional deconvolution.
 */
public class MitivDeconvolution extends EzPlug implements Block, EzStoppable {


    /***************************************************/
    /**         Plugin interface variables            **/
    /***************************************************/
    private EzVarBoolean expertMode;      // Tick box to expose expert parameters
    private EzTabs tabbedPane;

    /** data tab: **/
    private EzPanel  dataPanel;
    private EzVarSequence image;
    private EzVarChannel channel;
    private EzVarSequence psf; // Point Spread Function

    private EzVarInteger    paddingSizeXY, paddingSizeZ;
    private EzVarText imageSize, outputSize;
    /** weighting tab: **/
    private EzPanel   weigthPanel;
    private EzVarText weightsMethod;
    private final String[] weightOptions = new String[]{"None","Inverse covariance map","Variance map","Computed variance"};
    private EzVarDouble  gain, noise;
    private EzVarSequence weights, deadPixel;
    private EzButton showWeight;


    /** deconvolution tab: **/
    private EzPanel    deconvPanel;
    private EzVarDouble logmu, mu, epsilon;
    private EzVarSequence  restart;
    private EzVarBoolean  positivity;
    private EzVarBoolean  singlePrecision;
    private EzButton startDec, stopDec;
    private EzVarInteger    nbIterDeconv;

    /** headless mode: **/
    private EzVarSequence   outputHeadlessImage, outputHeadlessPSF;

    private EzLabel docDec;


    /****************************************************/
    /**                 VARIABLES                      **/
    /****************************************************/

    boolean debug =true;

    private int sizeX=512, sizeY=512, sizeZ=128; // Input sequence size
    protected int psfSizeX=1,psfSizeY=1,psfSizeZ=1;
    private Shape imageShape, psfShape;
    private  int Nxy=512, Nz=128;            // Output (padded sequence size)
    private Shape outputShape;

    Sequence cursequence=null;
    EdgePreservingDeconvolution solver ;
    private boolean run;
    private Vector result;

    /*********************************/
    /**      Initialization         **/
    /*********************************/
    @Override
    protected void initialize() {
        if (!isHeadLess()) {
            getUI().setParametersIOVisible(false);
            getUI().setActionPanelVisible(false);
        }


        tabbedPane = new EzTabs("BlindTabs", TabPlacement.TOP);

        /****************************************************/
        /**                    IMAGE TAB                   **/
        /****************************************************/
        //Creation of the inside of IMAGE TAB
        expertMode = new EzVarBoolean("Expert mode", false);

        dataPanel = new EzPanel("Step 1: Data"); //Border layout to be sure that the images are stacked to the up
        EzPanel imagePan = new EzPanel("FILEPanel");
        image = new EzVarSequence("Sequence:");
        channel = new EzVarChannel("Canal:", image.getVariable(), false);
        psf = new EzVarSequence("PSF:");
        imageSize = new EzVarText("Image size:");
        outputSize = new EzVarText("Output size:");
        paddingSizeXY = new EzVarInteger("padding xy:",0, Integer.MAX_VALUE,1);
        paddingSizeZ = new EzVarInteger("padding z :",0, Integer.MAX_VALUE,1);

        image.setNoSequenceSelection();
        psf.setNoSequenceSelection();

        expertMode.addVarChangeListener(new EzVarListener<Boolean>() {
            @Override
            public void variableChanged(EzVar<Boolean> source, Boolean newValue) {
                paddingSizeXY.setVisible(newValue);
                paddingSizeZ.setVisible(newValue);
                weigthPanel.setVisible(newValue);
                epsilon.setVisible(newValue);            }
        });

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



        image.addVarChangeListener(new EzVarListener<Sequence>() {
            @Override
            public void variableChanged(EzVar<Sequence> source,
                    Sequence newValue) {
                if (debug) {
                    System.out.println("Seq ch..."+image.getValue());
                }

                // getting metadata and computing sizes
                Sequence seq = image.getValue();
                if (seq != null || (seq != null && seq.isEmpty())) {
                    sizeX =  newValue.getSizeX();
                    sizeY = newValue.getSizeY();
                    sizeZ = newValue.getSizeZ();
                    updatePSFSize();
                    updateImageSize();

                    imageShape = new Shape(sizeX, sizeY, sizeZ);


                    if (debug) {
                        // show(IcyBufferedImageUtils.imageToArray(seq, imageShape,0),"Image map");
                        System.out.println("Seq changed:" + sizeX + "  "+ Nxy);
                    }
                    // setting restart value to the current sequence
                    restart.setValue(newValue);
                }
            }

        });

        psf.addVarChangeListener(new EzVarListener<Sequence>() {
            @Override
            public void variableChanged(EzVar<Sequence> source,
                    Sequence newValue) {
                if (debug) {
                    System.out.println("PSF changed"+psf.getValue());
                }
                // getting metadata and computing sizes
                //  Sequence seq = psf.getValue();
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









        // Note The listener of PSF is after BDEC tab

        /****************************************************/
        /**                WEIGHTING TAB                   **/
        /****************************************************/
        //Creation of the inside of WEIGHTING TAB
        weigthPanel = new EzPanel("Step 1b: Weights"); //Border layout to be sure that the images are stacked to the up
        EzPanel varianceTab = new EzPanel("VarianceTab");
        weightsMethod = new EzVarText(      "Weighting:", weightOptions, false);
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


        /****************************************************/
        /**                    DECONV TAB                  **/
        /****************************************************/
        //Creation of the inside of DECONVOLUTION TAB
        deconvPanel = new EzPanel("Step 2: Deconvolution"); //Border layout to be sure that the images are stacked to the up
        EzPanel deconvTab = new EzPanel("DeconvolutionTab");
        mu = new EzVarDouble("Regularization level:",1E-5,0.,Double.MAX_VALUE,0.01);
        logmu = new EzVarDouble("Log10 of the Regularization level:",-5,-Double.MAX_VALUE,Double.MAX_VALUE,1);
        epsilon = new EzVarDouble("Threshold level:",1E-2,0.,Double.MAX_VALUE,0.01);
        nbIterDeconv = new EzVarInteger("Number of iterations: ",10,0,Integer.MAX_VALUE ,1);
        positivity = new EzVarBoolean("Enforce nonnegativity:", true);
        // crop = new EzVarBoolean("Crop output to match input:", false);
        restart = new EzVarSequence("Starting point:");
        restart.setNoSequenceSelection();
        docDec = new EzLabel("Launch a deconvolution input PSF", Color.red);
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

        EzVarListener<Double> logmuActionListener = new EzVarListener<Double>() {
            @Override
            public void variableChanged(EzVar<Double> source, Double newValue) {
                mu.setValue(Math.pow(10, logmu.getValue()));
            }
        };
        logmu.addVarChangeListener(logmuActionListener);


        EzGroup groupStop1 = new EzGroup("Emergency STOP", stopDec);



        /**** IMAGE ****/
        imagePan.add(image);
        imagePan.add(channel);
        imagePan.add(imageSize);
        imagePan.add(psf);
        imagePan.add(paddingSizeXY);
        imagePan.add(paddingSizeZ);
        imagePan.add(outputSize);

        dataPanel.add(imagePan);
        tabbedPane.add(dataPanel);

        /**** Variance ****/
        varianceTab.add(weightsMethod);
        varianceTab.add(weights);
        varianceTab.add(gain);
        varianceTab.add(noise);
        varianceTab.add(deadPixel);
        varianceTab.add(showWeight);
        //Creation of VARIANCE TAB
        weigthPanel.add(varianceTab);
        tabbedPane.add(weigthPanel);

        /**** Deconv ****/

        deconvTab.add(logmu);
        deconvTab.add(mu);
        deconvTab.add(epsilon);
        deconvTab.add(nbIterDeconv);
        deconvTab.add(positivity);
        deconvTab.add(restart);
        deconvTab.add(docDec);
        deconvTab.add(startDec);
        deconvTab.add(groupStop1);
        //Creation of DECONVOLUTION TAB
        deconvPanel.add(deconvTab);
        tabbedPane.add(deconvPanel);

        addEzComponent(expertMode);
        addEzComponent(tabbedPane);
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
        //Here we will update the sequence
        //        if (sequence != null && sequence.isEmpty()){
        //            removeSequence(sequence);
        //        }

        if (sequence == null )  {

            sequence = new Sequence();
            addSequence(sequence);
        }

        if( sequence.getViewers() ==null){
            addSequence(sequence);
        }
        sequence.beginUpdate();

        Shape shape = arr.getShape();
        int rank = shape.rank();
        //      int type = arr.getType();
        int nx, ny, nz, nt;
        switch (rank) {
        case 2:
            nx  = shape.dimension(0);
            ny  = shape.dimension(1);
            sequence.setImage(0,0, new IcyBufferedImage(nx, ny, arr.flatten()));
            /*    switch (type) {
                    case Traits.DOUBLE:
                        sequence.setImage(0,0, new IcyBufferedImage(nx, ny, arr.toDouble().flatten()));
                        break;
                    case Traits.FLOAT:
                        sequence.setImage(0,0, new IcyBufferedImage(nx, ny, arr.toFloat().flatten()));
                        break;
                    default:
                        throwError("Show: only Double or Float array");
                        break;
                }*/
            break;
        case 3:
            nx  = shape.dimension(0);
            ny  = shape.dimension(1);
            nz =  shape.dimension(2);


            for (int j = 0; j < nz; j++) {
                sequence.setImage(0,j, new IcyBufferedImage(nx, ny,((Array3D)arr).slice(j).flatten() ));
            }/*
                switch (type) {
                    case Traits.DOUBLE:
                        for (int j = 0; j < nz; j++) {
                            sequence.setImage(0,j, new IcyBufferedImage(nx, ny, ((Double3D) arr).slice(j).flatten()));
                        }
                        break;
                    case Traits.FLOAT:
                        for (int j = 0; j < nz; j++) {
                            sequence.setImage(0,j, new IcyBufferedImage(nx, ny, ((Float3D) arr).slice(j).flatten()));
                        }
                        break;
                    default:
                        throwError("Show: only Double or Float array");
                        break;
                }*/

            break;

        case 4:
            nx  = shape.dimension(0);
            ny  = shape.dimension(1);
            nz =  shape.dimension(2);
            nt =  shape.dimension(3);

            for (int k = 0; k < nt; k++) {
                for (int j = 0; j < nz; j++) {
                    sequence.setImage(k,j, new IcyBufferedImage(nx, ny,((Array3D)arr).slice(k).slice(j).flatten() ));
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

        int numCanal = channel.getValue();

        DoubleArray imgArray, psfArray, wgtArray;

        imgArray = (DoubleArray) IcyBufferedImageUtils.imageToArray(imgSeq, imageShape, numCanal);
        psfArray = (DoubleArray) IcyBufferedImageUtils.imageToArray(psfSeq, psfShape, 0);
        wgtArray = createWeights(imgArray).toDouble();

        cursequence = new Sequence("Current Iterate");
        cursequence.copyMetaDataFrom(imgSeq, false);

        addSequence(cursequence);



        solver.setObjectShape(outputShape);
        solver.setPSF(psfArray);
        solver.setData(imgArray);
        solver.setWeight(wgtArray);
        solver.setEdgeThreshold(epsilon.getValue());
        solver.setRegularizationLevel(mu.getValue());
        //    solver.setForceSinglePrecision(singlePrecision.getValue());
        solver.setSaveBest(true);
        if (positivity.getValue()){
            solver.setLowerBound(0.);
        }
        solver.setMaximumIterations(nbIterDeconv.getValue());
        solver.setMaximumEvaluations(2*nbIterDeconv.getValue());
        System.out.println("Launch it:"+nbIterDeconv.getValue());
        run = true;
        OptimTask task = solver.start();
        show(solver.getPSF(),"PSF");
        show(solver.getData(),"Data");

        while (run) {
            task = solver.getTask();
            System.out.println("  it "+solver.getIterations());
            if (task == OptimTask.ERROR) {
                System.err.format("Error: %s\n", solver.getReason());
                break;
            }
            if (task == OptimTask.NEW_X || task == OptimTask.FINAL_X) {
                show(solver.getSolution(),cursequence,"Current mu="+solver.getRegularizationLevel() +"it:"+solver.getIterations());
                if (task == OptimTask.FINAL_X) {
                    break;
                }
            }
            if (task == OptimTask.WARNING) {
                break;
            }
            solver.iterate();
        }
        //   Vector vect;

        //   DoubleArray dblArr = DoubleArray.createFrom( solver.getBestSolution());
        // show(solver.getBestSolution(),"Result ");

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
        inputMap.add("psf", psf.getVariable());
        inputMap.add("mu", mu.getVariable());
        inputMap.add("epsilon", epsilon.getVariable());
        inputMap.add("maxIter", nbIterDeconv.getVariable());
        inputMap.add("Padx", paddingSizeXY.getVariable());
        inputMap.add("PadZ", paddingSizeZ.getVariable());
    }

    //The output variable for the protocol
    @Override
    public void declareOutput(VarList outputMap) {
        //  outputMap.add("output", outputMap.getVariable());
    }

    @Override
    public void clean() {
        // TODO Auto-generated method stub

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