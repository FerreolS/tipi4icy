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


package plugins.ferreol.demics;

import static plugins.mitiv.io.Icy2TiPi.arrayToSequence;
import static plugins.mitiv.io.Icy2TiPi.sequenceToArray;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.SwingUtilities;

import icy.file.Loader;
import icy.gui.frame.progress.AnnounceFrame;
import icy.image.colormap.IceColorMap;
import icy.main.Icy;
import icy.sequence.Sequence;
import microTiPi.epifluorescence.WideFieldModel;
import microTiPi.microscopy.PSF_Estimation;
import mitiv.array.ArrayUtils;
import mitiv.array.Double2D;
import mitiv.array.DoubleArray;
import mitiv.base.Shape;
import mitiv.cost.EdgePreservingDeconvolution;
import mitiv.optim.OptimTask;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzButton;
import plugins.adufour.ezplug.EzGroup;
import plugins.adufour.ezplug.EzLabel;
import plugins.adufour.ezplug.EzPanel;
import plugins.adufour.ezplug.EzStoppable;
import plugins.adufour.ezplug.EzTabs;
import plugins.adufour.ezplug.EzTabs.TabPlacement;
import plugins.adufour.ezplug.EzVar;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarChannel;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarDoubleArrayNative;
import plugins.adufour.ezplug.EzVarFile;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarListener;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.ezplug.EzVarText;

/**
 * This class implements  EpiDEMIC, an Icy plugin for 3D blind deconvolution in epifluorescence (wide field) fluorescence microscopy.
 *
 * @author Ferréol Soulez & Jonathan Leger
 *
 */
public class EpiDEMIC extends DEMICSPlug implements  EzStoppable, Block {


    /***************************************************/
    /**         Plugin interface variables            **/
    /***************************************************/
    private EzTabs tabbedPane;              // The interface is composed of several tabs

    /** data tab: **/
    private EzPanel         dataPanel;      // data tab
    private EzButton        saveMetaData, showPSF;

    /** weighting tab: **/
    private EzGroup         ezWeightingGroup, ezPadGroup;


    /** deconvolution tab: **/
    private EzPanel         deconvPanel;
    private EzVarDouble     epsilon; // deconvolution hyper parameters; mu = 10^(logmu)

    private EzGroup         ezDeconvolutionGroup;
    protected  int          Nxy=128; // Output (padded sequence size)


    /** blind deconvolution tab: **/
    private EzPanel         bdecPanel;
    private EzVarText       nbAlphaCoef;    // number of mode to describe the phase
    private EzVarText       nbBetaCoef;     // number of mode to describe the modulus
    private final String[]  nAlphaOptions = new String[]{"0","1","3","7","8","12","18","19","25","33","34","42","52","53","63","75","76","88","102","103"};
    private final String[]  nBetaOptions = new String[]{"0","3","4","6","10","11","15","21","22","28","36","37","45","55","56","66","78","79","91","105","106"};
    private final String[]  nAlphaOptionsR = new String[]{"0","1","2","3","4","5","6","7","8","9"};
    private final String[]  nBetaOptionsR = new String[]{"0","1","2","3","4","5","6","7","8","9"};
    private EzVarBoolean    radial;         // use only radial mode (constraint the PSF to be radially symmetric)
    private EzButton        showPSF2, showModulus, showPhase;
    private EzButton        startBlind, stopBlind, showFullObject2, resetPSF;
    private EzButton        saveParam, loadParam;
    private EzLabel         docDec;
    private EzVarInteger    totalNbOfBlindDecLoop;// number of outer loop
    private EzVarInteger    maxIterDefocus; // number of iteration for defocus parameters
    private EzVarInteger    maxIterPhase;   // number of iteration for phase parameters
    private EzVarInteger    maxIterModulus; // number of iteration for modulus parameters
    private EzLabel         docBlind;
    private EzGroup         visuPSF;
    private EzGroup         ezBlindDeconvolutionGroup;

    // PSF parameters
    private EzVarDoubleArrayNative pupilShift;  // estimated shift of the pupil relatively to the pupil axis
    private EzVarDoubleArrayNative phaseCoefs, modulusCoefs; // estimated coefs of phase and modulus

    /** headless mode: **/
    private EzVarSequence   outputHeadlessPSF=null;

    // Global solvers
    private PSF_Estimation psfEstimation;
    private EdgePreservingDeconvolution solver;



    // Main arrays for the psf estimation
    private boolean runBdec;
    private WideFieldModel pupil=null;
    private boolean guessModulus;
    private boolean guessPhase;


    /*********************************/
    /**            DEBUG            **/
    /*********************************/
    private boolean debug = false;      // Show psf steps
    private boolean verbose = false;    // Show some values, need debug to true
    private EzPanel  debuggingPanel;
    private EzVarText resultCostPrior, resultDefocus, resultPhase, resultModulus;


    private String psfPath=null;





    @Override //FIXME should be in another class
    public void clean() {
    }





    /**
     * Update the size of the deconvolved image according the size of the input and the padding rounded to the next best fft size
     *
     */
    @Override
    protected void updatePaddedSize() {
        super.updatePaddedSize();
        Nxy = Math.max(Nx, Ny);
        psfShape = new Shape(Nxy, Nxy, Nz);
        outputShape = new Shape(Nxy, Nxy, Nz);

    }


    /**
     *  set default values of the plugin
     */
    @Override
    protected void setDefaultValue() {
        super.setDefaultValue();
        radial.setValue(false);
        paddingSizeXY.setValue(30);

        if(debug){
            resultCostPrior.setValue(   "No results yet");
            resultDefocus.setValue(     "No results yet");
            resultModulus.setValue(     "No results yet");
            resultPhase.setValue(       "No results yet");
        }

        if (!isHeadLess()) {
            if(debug){
                resultCostPrior.setEnabled(false);
                resultDefocus.setEnabled(false);
                resultModulus.setEnabled(false);
                resultPhase.setEnabled(false);
            }
        }
    }


    /* (non-Javadoc)
     *      Initialization  of the plugin
     * @see plugins.adufour.ezplug.EzPlug#initialize()
     */
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

        dataPanel = new EzPanel("Step 1: Data");
        data = new EzVarSequence("Sequence:");
        channel = new EzVarChannel("Channel:", data.getVariable(), false);
        dataSize = new EzVarText("Image size:");
        outputSize = new EzVarText("Output size:");
        paddingSizeXY = new EzVarInteger("padding xy:",0, Integer.MAX_VALUE,1);
        paddingSizeZ = new EzVarInteger("padding z :",0, Integer.MAX_VALUE,1);
        restart = new EzVarSequence("Starting point:");
        restart.setNoSequenceSelection();
        channelRestart = new EzVarChannel("Initialization channel :", restart.getVariable(), false);
        saveMetaData = new EzButton("Save metadata", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if (debug) {
                    System.out.println("Saving metadata");
                }
                updateMetaData();
            }
        });



        EzVarListener<Integer> zeroPadActionListener = new EzVarListener<Integer>() {
            @Override
            public void variableChanged(EzVar<Integer> source, Integer newValue) {
                updatePaddedSize();
                updateImageSize();
                updateOutputSize();
                pupil=null;
            }
        };
        paddingSizeXY.addVarChangeListener(zeroPadActionListener);
        paddingSizeZ.addVarChangeListener(zeroPadActionListener);



        data.addVarChangeListener(new EzVarListener<Sequence>() {
            @Override
            public void variableChanged(EzVar<Sequence> source,
                    Sequence newValue) {
                if (debug) {
                    System.out.println("Seq ch..."+newValue);
                }
                dataChanged() ;
                if(dataSeq!=null){
                    meta = getMetaData(dataSeq);
                    dxy_nm.setValue(    meta.dxy);
                    dz_nm.setValue(     meta.dz);
                    scale.setValue(new double[]{1.0 ,1.0, dxy_nm.getValue()/ dz_nm.getValue() } );
                    na.setValue(     meta.na);
                    lambda.setValue( meta.lambda);
                    ni.setValue(     meta.ni);
                    if (debug) {
                        System.out.println("Seq changed:" + sizeX + "  "+ Nxy);
                    }
                    // setting restart value to the current sequence
                    restart.setValue(newValue);
                }

            }

        });

        channel.addVarChangeListener(new EzVarListener<Integer>() {
            @Override
            public void variableChanged(EzVar<Integer> source, Integer newValue) {
                channelRestart.setValue(newValue);
            }
        });
        restart.addVarChangeListener(new EzVarListener<Sequence>() {
            @Override
            public void variableChanged(EzVar<Sequence> source,
                    Sequence newValue) {
                newValue = restart.getValue();
                if(debug){
                    if (newValue != null || (newValue != null && newValue.isEmpty())) {
                        System.out.println("restart changed:"+newValue.getName());
                    }
                }
            }
        });
        EzVarListener<Double> metaActionListener = new EzVarListener<Double>() {
            @Override
            public void variableChanged(EzVar<Double> source, Double newValue) {
                scale.setValue(new double[]{1.0 ,1.0,  dz_nm.getValue()/dxy_nm.getValue() } );
                //pupilShift.setValue(new double[] { 0., 0.});
                pupil=null;
            };
        };



        dxy_nm = new EzVarDouble("dxy(nm):",64.5,0., Double.MAX_VALUE,1.);
        dxy_nm.addVarChangeListener(metaActionListener);
        dz_nm = new EzVarDouble("dz(nm):",128.,0., Double.MAX_VALUE,1.);
        dz_nm.addVarChangeListener(metaActionListener);

        na = new EzVarDouble("NA:",1.4,0.,Double.MAX_VALUE,0.05);
        na.addVarChangeListener(metaActionListener);
        ni = new EzVarDouble("ni:",1.518,1.,2.,0.01);
        ni.addVarChangeListener(metaActionListener);
        lambda = new EzVarDouble( "\u03BB(nm):",540.,10.,15000.,5);
        lambda.addVarChangeListener(metaActionListener);

        showPSF = new EzButton("Show PSF", new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                // Show the initial PSF
                psfClicked();
                if (debug) {
                    System.out.println("Showing PSF Button");
                }
            }
        });


        /****************************************************/
        /**                WEIGHTING Group                   **/
        /****************************************************/
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
                } else if (weightsMethod.getValue() == weightOptions[1] || weightsMethod.getValue() == weightOptions[2]) {  //Personalized map or Variance map
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
        showWeight = new EzButton("Show weight map", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                showWeightClicked();
                if (debug) {
                    System.out.println("Weight compute");
                }
            }
        });

        ezPadGroup = new EzGroup("Padding",paddingSizeXY, paddingSizeZ);
        ezPadGroup.setFoldedState(true);
        ezWeightingGroup = new EzGroup("Weighting",weightsMethod,weights,gain,noise,deadPixel,showWeight);
        ezWeightingGroup.setFoldedState(true);

        loadFile = new EzVarFile("Load parameters from", "","*.xml");
        loadParam = new EzButton("Load parameters", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(loadFile.getValue()!=null){
                    loadParamClicked();
                }
                if (debug) {
                    System.out.println("Load parameters");
                }
            }
        });

        /****************************************************/
        /**                    DECONV TAB                  **/
        /****************************************************/
        deconvPanel = new EzPanel("Step 2: Deconvolution");
        mu = new EzVarDouble("Regularization level:",1E-5,0.,Double.MAX_VALUE,0.01);
        logmu = new EzVarDouble("Log10 of the Regularization level:",-5,-Double.MAX_VALUE,Double.MAX_VALUE,1);
        epsilon = new EzVarDouble("Threshold level:",1E-2,0.,Double.MAX_VALUE,1.0);
        scale = new EzVarDoubleArrayNative("Aspect ratio of a voxel", new double[][] { {1.0 ,1.0, 1.0} }, true);
        nbIterDeconv = new EzVarInteger("Number of iterations: ",10,0,Integer.MAX_VALUE ,1);
        positivity = new EzVarBoolean("Enforce nonnegativity:", true);
        singlePrecision = new EzVarBoolean("Compute in single precision:", false);
        singlePrecision.addVarChangeListener(new EzVarListener<Boolean>() {
            @Override
            public void variableChanged(EzVar<Boolean> source, Boolean newValue) {
                pupil=null;
                psfEstimation=null;
                cursequence =null;

            }
        });
        docDec = new EzLabel("Launch a deconvolution using the current PSF", Color.red);
        startDec = new EzButton("Start Deconvolution", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                Thread workerThread = new Thread() {
                    @Override
                    public void run() {
                        launch(true);

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

        showFullObject = new EzButton("Show the full (padded) object", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(debug){
                    System.out.println("showFull");
                }
                Sequence fSeq;
                fSeq = new Sequence("Deconvolved image");
                fSeq.copyMetaDataFrom(data.getValue(), false);
                if(objArray != null){
                    if(solver != null){
                        show(objArray,fSeq,"Deconvolved "+ data.getValue().getName() + " with padding. mu: " +solver.getRegularizationLevel() );
                    }
                    else{
                        show(objArray,fSeq,"Deconvolved "+ data.getValue().getName() + " with padding. mu " + mu.getValue());
                    }
                }else if(solver != null){
                    show(solver.getSolution(),fSeq,"Deconvolved "+ data.getValue().getName() + "with padding. mu="+solver.getRegularizationLevel() );
                }else {
                    show(ArrayUtils.extract(dataArray, outputShape),fSeq,"Deconvolved "+ data.getValue().getName() + "with padding. mu="+ mu.getValue() );
                }
            }
        });

        EzGroup groupStop1 = new EzGroup("Emergency STOP", stopDec);

        ezDeconvolutionGroup = new EzGroup("Expert  parameters",epsilon,scale,positivity,singlePrecision, showFullObject);
        ezDeconvolutionGroup.setFoldedState(true);

        /****************************************************/
        /**                      BDEC TAB                  **/
        /****************************************************/
        //Saving variables
        pupilShift = new EzVarDoubleArrayNative("pupilShift", new double[][] { { 0.0, 0.0} }, false);
        pupilShift.setVisible(false);
        phaseCoefs = new EzVarDoubleArrayNative("phase coefs",null , false);
        phaseCoefs.setVisible(false);
        modulusCoefs = new EzVarDoubleArrayNative("modulusCoefs", new double[][] { {1.0 } },0, false);
        modulusCoefs.setVisible(false);

        bdecPanel = new EzPanel("Step 3: Blind dec.");
        nbAlphaCoef = new EzVarText(            "Number of phase coefs N\u03B1:", nAlphaOptions, 3,false );
        nbBetaCoef = new EzVarText(             "Number of modulus coefs N\u03B2:", nBetaOptions,0, false);
        radial = new EzVarBoolean(              "Radially symmetric PSF", false);
        maxIterDefocus = new EzVarInteger(      "Max. nb. of iterations for defocus:",20,0,Integer.MAX_VALUE ,1);
        maxIterPhase = new EzVarInteger(        "Max. nb. of iterations for phase:",20,0,Integer.MAX_VALUE ,1);
        maxIterModulus = new EzVarInteger(      "Max. nb. of iterations for modulus:",0,0,Integer.MAX_VALUE ,1);
        totalNbOfBlindDecLoop = new EzVarInteger(  "Number of loops:",2,0,Integer.MAX_VALUE ,1);
        showIteration = new EzVarBoolean("Show intermediate results:", true);

        if (isHeadLess()){
            showIteration.setValue(false);
        }


        resetPSF = new EzButton(            "Reset PSF", new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                pupil=null;
            }
        });

        nbAlphaCoef.addVarChangeListener(new EzVarListener<String>() {
            @Override
            public void variableChanged(EzVar<String> source, String newValue) {
                pupil=null;
            }
        });


        nbBetaCoef.addVarChangeListener(new EzVarListener<String>() {
            @Override
            public void variableChanged(EzVar<String> source, String newValue) {
                pupil=null;
            }
        });

        radial.addVarChangeListener(new EzVarListener<Boolean>() {
            @Override
            public void variableChanged(EzVar<Boolean> source, Boolean newValue) {
                if(newValue){
                    nbAlphaCoef.setDefaultValues( nAlphaOptionsR,3, false );
                    nbBetaCoef.setDefaultValues( nBetaOptionsR,1, false );
                    pupil=null;
                }else{
                    nbAlphaCoef.setDefaultValues( nAlphaOptions,7, false );
                    nbBetaCoef.setDefaultValues(  nBetaOptions,1, false );
                    pupil=null;
                }

            }
        });

        showPSF2 = new EzButton(        "Show PSF", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                psfClicked();
                if (debug) {
                    System.out.println("First PSF compute");
                }
            }
        });

        showPhase = new EzButton(       "Phase of the pupil function", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                phaseClicked();
                if (debug) {
                    System.out.println("Show phase");
                }
            }
        });

        showModulus = new EzButton(     "Modulus of the pupil function", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                modulusClicked();
                if (debug) {
                    System.out.println("Show modulus");
                }
            }
        });
        docBlind = new EzLabel("Joint estimation of the image and PSF", Color.RED);
        startBlind = new EzButton("Blind deconvolution", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                Thread workerThread = new Thread() {
                    @Override
                    public void run() {
                        long startTime = System.currentTimeMillis();
                        if (data.getValue()!=null)
                            launch(false);


                        long stopTime = System.currentTimeMillis();
                        long elapsedTime = stopTime - startTime;
                        System.out.println("time: "+elapsedTime);
                    }
                };
                workerThread.start();
            }
        });

        stopBlind = new EzButton("STOP Computation", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                stopExecution();
            }
        });

        EzGroup groupStop2 = new EzGroup("Emergency STOP", stopBlind);

        showFullObject2 = new EzButton("Show the full (padded) object", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(debug){
                    System.out.println("showFull");
                }
                Sequence fSeq;
                fSeq = new Sequence("Deconvolved image");
                fSeq.copyMetaDataFrom(data.getValue(), false);
                if(objArray != null){
                    if(solver != null){
                        show(objArray,fSeq,"Deconvolved "+ data.getValue().getName() + " with padding. mu: " +solver.getRegularizationLevel() );
                    }
                    else{
                        show(objArray,fSeq,"Deconvolved "+ data.getValue().getName() + " with padding. mu " + mu.getValue());
                    }
                }else if(solver != null){
                    show(solver.getSolution(),fSeq,"Deconvolved "+ data.getValue().getName() + "with padding. mu="+solver.getRegularizationLevel() );
                }else {
                    show(ArrayUtils.extract(dataArray, outputShape),fSeq,"Deconvolved "+ data.getValue().getName() + "with padding. mu="+ mu.getValue() );
                }
            }
        });



        visuPSF = new EzGroup("Visualization", showFullObject2,showIteration,showPSF2, showPhase, showModulus);

        saveFile = new EzVarFile("Save parameters in", "");
        saveParam = new EzButton("Save parameters", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                saveParamClicked();
                if (debug) {
                    System.out.println("Save parameters");
                }
            }
        });

        ezBlindDeconvolutionGroup = new EzGroup("Expert  parameters",nbAlphaCoef,nbBetaCoef,
                radial, maxIterDefocus, maxIterPhase, maxIterModulus);
        ezBlindDeconvolutionGroup.setFoldedState(true);

        /****************************************************/
        /**                    RESULT TAB                  **/
        /****************************************************/
        debuggingPanel = new EzPanel("Results");
        resultCostPrior = new EzVarText(    "Costs");
        resultDefocus = new EzVarText(      "Defocus");
        resultModulus = new EzVarText(      "Modulus");
        resultPhase = new EzVarText(        "Phase");



        /****************************************************/
        /**                      ToolTips                  **/
        /****************************************************/
        data.setToolTipText(ToolTipText.sequenceImage);
        weights.setToolTipText(ToolTipText.sequenceWeigth);
        weightsMethod.setToolTipText(ToolTipText.sequenceWeigth);
        deadPixel.setToolTipText(ToolTipText.sequencePixel);

        channel.setToolTipText(ToolTipText.textCanal);

        dxy_nm.setToolTipText(ToolTipText.doubleDxy);
        dz_nm.setToolTipText(ToolTipText.doubleDz);
        na.setToolTipText(ToolTipText.doubleNa);
        lambda.setToolTipText(ToolTipText.doubleLambda);
        ni.setToolTipText(ToolTipText.doubleNi);
        gain.setToolTipText(ToolTipText.doubleGain);
        noise.setToolTipText(ToolTipText.doubleNoise);
        mu.setToolTipText(ToolTipText.doubleMu);
        epsilon.setToolTipText(ToolTipText.doubleEpsilon);
        nbIterDeconv.setToolTipText(ToolTipText.doubleMaxIter);
        paddingSizeXY.setToolTipText(ToolTipText.doublePadding);
        paddingSizeZ.setToolTipText(ToolTipText.doublePadding);

        nbAlphaCoef.setToolTipText(ToolTipText.doubleNalpha);
        nbBetaCoef.setToolTipText(ToolTipText.doubleNbeta);
        totalNbOfBlindDecLoop.setToolTipText(ToolTipText.doubleBDecTotalIteration);

        restart.setToolTipText(ToolTipText.booleanRestart);
        positivity.setToolTipText(ToolTipText.booleanPositivity);
        showFullObject.setToolTipText(ToolTipText.booleanCrop);
        showFullObject2.setToolTipText(ToolTipText.booleanCrop);
        debuggingPanel.setToolTipText(ToolTipText.textOutput);

        /******** Adding ********/
        /*   paddingSizeXY.setVisible(false);
        paddingSizeZ.setVisible(false);
        epsilon.setVisible(false);
        scale.setVisible(false);
        singlePrecision.setVisible(false);

        nbAlphaCoef.setVisible(false);
        nbBetaCoef.setVisible(false);
        radial.setVisible(false);
        maxIterDefocus.setVisible(false);
        maxIterPhase.setVisible(false);
        maxIterModulus.setVisible(false);*/

        /**** IMAGE ****/
        dataPanel.add(data);
        dataPanel.add(channel);
        dataPanel.add(dataSize);
        dataPanel.add(ezPadGroup);
        dataPanel.add(outputSize);
        dataPanel.add(dxy_nm);
        dataPanel.add(dz_nm);


        dataPanel.add(na);
        dataPanel.add(ni);
        dataPanel.add(lambda);

        dataPanel.add(ezWeightingGroup);

        dataPanel.add(saveMetaData);
        dataPanel.add(showPSF);
        dataPanel.add(loadFile);
        dataPanel.add(loadParam);
        tabbedPane.add(dataPanel);


        /**** Deconv ****/

        deconvPanel.add(logmu);
        deconvPanel.add(mu);
        deconvPanel.add(nbIterDeconv);
        deconvPanel.add(restart);
        deconvPanel.add(channelRestart);
        deconvPanel.add(ezDeconvolutionGroup);

        deconvPanel.add(docDec);
        deconvPanel.add(startDec);
        deconvPanel.add(groupStop1);
        tabbedPane.add(deconvPanel);

        /**** BDec ****/
        bdecPanel.add(resetPSF);
        bdecPanel.add(totalNbOfBlindDecLoop);
        bdecPanel.add(ezBlindDeconvolutionGroup);
        bdecPanel.add(docBlind);
        bdecPanel.add(startBlind);
        bdecPanel.add(visuPSF);
        bdecPanel.add(groupStop2);
        bdecPanel.add(saveFile);
        bdecPanel.add(saveParam);
        tabbedPane.add(bdecPanel);
        if(debug){
            /**** Result ****/
            debuggingPanel.add(resultCostPrior);
            debuggingPanel.add(resultDefocus);
            debuggingPanel.add(resultModulus);
            debuggingPanel.add(resultPhase);
            tabbedPane.add(debuggingPanel);
        }

        addEzComponent(pupilShift);
        addEzComponent(phaseCoefs);
        addEzComponent(modulusCoefs);
        addEzComponent(tabbedPane);
        // Must be added to global panel first
        visuPSF.setFoldedState(true);


        setDefaultValue();

        updatePaddedSize();
        updateOutputSize();
        updateImageSize();

        if (isHeadLess()) {
            outputHeadlessImage = new EzVarSequence("Output Image");
            outputHeadlessPSF = new EzVarSequence("Output PSF");
            outputHeadlessWght = new EzVarSequence("Computed weight");
        }

    }

    /**
     * Load the parameter file and perform parameter update
     */
    private void loadParamClicked() {
        File saveName = saveFile.getValue();
        this.loadParameters(loadFile.getValue());
        saveFile.setValue(saveName);    // FIX saving file name erasing during load

        buildpupil();
        pupil.setPupilAxis(pupilShift.getValue());
        pupil.setNi(ni.getValue());
        pupil.setModulus(modulusCoefs.getValue());
        pupil.setPhase(phaseCoefs.getValue());
        if (dataSeq != null) {
            sizeX = dataSeq.getSizeX();
            sizeY = dataSeq.getSizeY();
            sizeZ = dataSeq.getSizeZ();
            if (sizeZ == 1) {
                throwError("Input data must be 3D");
                return;
            }
            updatePaddedSize();
            updateOutputSize();
            updateImageSize();
            dataShape = new Shape(sizeX, sizeY, sizeZ);
        }

    }

    private void saveParamClicked() {
        if(pupil!=null){
            pupilShift.setValue( pupil.getPupilShift());
            if(pupil.getPhaseCoefs() !=null)
                phaseCoefs.setValue(pupil.getPhaseCoefs().getData());
            modulusCoefs.setValue(pupil.getModulusCoefs().getData());
        }
        if (debug) {
            System.out.println("--------------");
            System.out.println("defocus");
            System.out.println( pupil.getDefocus().toString() );
        }
        File pathName = saveFile.getValue();
        if(pathName!=null){
            if (!pathName.getName().endsWith(".xml")){
                pathName = new File(pathName.getAbsolutePath()+".xml");
            }

            this.saveParameters(pathName);
        }
    }


    @Override
    protected void execute() {

        if (isHeadLess()){

            if(  Icy.getCommandLinePluginArgs().length!=0){
                initialize();
                parseCmdLine();
            }
            showIteration.setValue(false);
        }
        long startTime = System.currentTimeMillis();
        launch(false);

        long stopTime = System.currentTimeMillis();
        long elapsedTime = stopTime - startTime;
        System.out.println("time: "+elapsedTime);
    }

    private void launch(boolean runDeconv) {
        try {
            startBlind.setText("Computing...");

            buildpupil();
            if (debug|| isHeadLess()) {
                System.out.println("-------------IMAGE-------------------");
                System.out.println("File: "+data.getValue());
                System.out.println("Canal: "+channel.getValue());
                System.out.println("image size: "+ dataSize.getValue());
                System.out.println("--------------PSF------------------");
                System.out.println("dxy: "+dxy_nm.getValue()*1E-9);
                System.out.println("dz: "+dz_nm.getValue()*1E-9);
                System.out.println("Nxy: "+Nxy);
                System.out.println("Nx: "+Nz);
                System.out.println("NA: "+na.getValue());
                System.out.println("\u03BB: "+lambda.getValue()*1E-9);
                System.out.println("ni: "+ni.getValue());
                System.out.println("--------------Variance------------------");
                System.out.println("Weights method: "+weightsMethod.getValue());
                System.out.println("Weights: "+weights.getValue());
                System.out.println("Gain: "+gain.getValue());
                System.out.println("Noise: "+noise.getValue());
                System.out.println("deadPix: "+deadPixel.getValue());
                System.out.println("--------------DECONV------------------");

                System.out.println("zeroPad xy: "+paddingSizeXY.getValue());
                System.out.println("zeroPad z: "+paddingSizeZ.getValue());
                System.out.println("nbIter: "+nbIterDeconv.getValue());
                System.out.println("Restart: "+restart.getValue());
                System.out.println("Positivity: "+positivity.getValue());
                System.out.println("--------------BDEC------------------");
                System.out.println("output size: "+ outputSize.getValue());
                System.out.println("nbIter: "+nbIterDeconv.getValue());
                System.out.println("Number of total iterations: "+totalNbOfBlindDecLoop.getValue());
                System.out.println("------------------------------------");
                System.out.println("");
            }
            run = true;
            runBdec = !runDeconv;

            preProcessing();

            /*---------------------------------------*/
            /*            OPTIMISATION               */
            /*---------------------------------------*/

            if (runBdec) {
                if ((pupil.getPhaseCoefs()==null)||( Integer.parseInt(nbAlphaCoef.getValue()) != pupil.getNPhase())){
                    int nbAlpha = Integer.parseInt(nbAlphaCoef.getValue());
                    if (nbAlpha==0){
                        guessPhase = false;
                    }else{
                        guessPhase = true;
                        pupil.setNPhase(nbAlpha);
                    }
                }

                if  ((pupil.getModulusCoefs()==null)||(Integer.parseInt(nbBetaCoef.getValue()) != pupil.getNModulus())){
                    int nbBeta = Integer.parseInt(nbBetaCoef.getValue());


                    if (nbBeta==0){
                        guessModulus = false;
                    }else{
                        guessModulus = true;
                        pupil.setNModulus(nbBeta);
                    }
                }


                psfEstimation = new PSF_Estimation(pupil);

                psfEstimation.setWeight(  ArrayUtils.pad(wgtArray,outputShape));
                psfEstimation.setData(ArrayUtils.pad(dataArray,outputShape));

                psfEstimation.enablePositivity(false);
                psfEstimation.setAbsoluteTolerance(0.0);

                for(int i = 0; i < totalNbOfBlindDecLoop.getValue(); i++) {
                    psfArray = ArrayUtils.roll(pupil.getPSF());
                    pupil.freePSF();
                    deconv();

                    psfEstimation.setObj(objArray);

                    /* Defocus estimation */
                    if (maxIterDefocus.getValue()>0){
                        if (debug && verbose) {
                            System.out.println("------------------");
                            System.out.println("Defocus estimation");
                            System.out.println("------------------");
                        }
                        psfEstimation.setRelativeTolerance(0.);
                        psfEstimation.setMaximumIterations(maxIterDefocus.getValue());
                        psfEstimation.fitPSF( WideFieldModel.DEFOCUS);
                    }

                    /* Phase estimation */
                    if((maxIterPhase.getValue()>0)& guessPhase){
                        if (debug && verbose) {
                            System.out.println("Phase estimation");
                            System.out.println("------------------");
                        }
                        psfEstimation.setMaximumIterations(maxIterPhase.getValue());
                        psfEstimation.fitPSF( WideFieldModel.PHASE);
                    }


                    /* Modulus estimation */
                    if((maxIterModulus.getValue() >0)&guessModulus){
                        if (debug && verbose) {
                            System.out.println("Modulus estimation");
                            System.out.println("------------------");
                        }
                        psfEstimation.setMaximumIterations(maxIterModulus.getValue());
                        psfEstimation.fitPSF( WideFieldModel.MODULUS);
                    }

                    //Emergency stop
                    if (!run) {
                        return;
                    }
                }
            } else {
                psfArray = ArrayUtils.roll( pupil.getPSF() );
                pupil.freePSF();
                preProcessing();
                deconv();
            }
            if(pupil!=null)
                pupil.freePSF();// TODO free more memory

            SwingUtilities.invokeLater(new Runnable() {
                @Override
                public void run() {
                    restart.setValue(cursequence);
                    channelRestart.setValue(0);
                    ni.setValue(pupil.getNi());
                    if (isHeadLess()) {
                        if(outputHeadlessImage==null){
                            outputHeadlessImage = new EzVarSequence("Output Image");
                        }
                        if(outputHeadlessPSF==null){
                            outputHeadlessPSF = new EzVarSequence("Output PSF");
                        }

                        if (outputHeadlessWght==null) {
                            outputHeadlessWght = new EzVarSequence("Computed weights");
                        }

                        Sequence psfSequence = null;
                        psfSequence =   arrayToSequence( ArrayUtils.roll(pupil.getPSF()), psfSequence);


                        outputHeadlessPSF.setValue(psfSequence);
                        outputHeadlessImage.setValue(cursequence);
                        outputHeadlessWght.setValue(arrayToSequence(wgtArray));

                        if(outputPath!=null){
                            saveSequence(cursequence, outputPath);
                        }

                        if(psfPath!=null){
                            saveSequence(psfSequence, psfPath);
                        }

                        if(saveFile.getValue()!=null){
                            saveParamClicked();
                        }
                    }
                }
            });
        } catch (IllegalArgumentException e) {
            new AnnounceFrame("Oops, Error: "+ e.getMessage());
            if (debug) {
                e.printStackTrace();
            }
        } finally {
            startBlind.setText("Guess PSF");
        }
    }


    /*****************************************/
    /** All the PSF buttons call are here   **/
    /*****************************************/

    private void buildpupil()
    {

        if(pupil==null)
        {
            pupil = new WideFieldModel(psfShape,Integer.parseInt(nbAlphaCoef.getValue()),
                    Integer.parseInt(nbBetaCoef.getValue()), na.getValue(),
                    lambda.getValue()*1E-9, ni.getValue(),
                    dxy_nm.getValue()*1E-9, dz_nm.getValue()*1E-9,
                    radial.getValue(), singlePrecision.getValue());
            pupil.setPupilAxis(pupilShift.getValue());
            pupil.setModulus(modulusCoefs.getValue());
            if(phaseCoefs.getValue()!=null)
                pupil.setPhase(phaseCoefs.getValue());
        }
    }

    private void psfClicked()
    {   // Compute the PSF in a Thread to prevent GUI freezing
        Thread workerThread = new Thread() {
            @Override
            public void run() {

                buildpupil();
                pupil.computePSF();

                SwingUtilities.invokeLater(new Runnable() {
                    @Override
                    public void run() {
                        Sequence psfSequence = new Sequence();
                        if(dataSeq!=null){
                            psfSequence.copyMetaDataFrom(dataSeq, false);
                        }
                        show(ArrayUtils.roll(pupil.getPSF()),psfSequence,"Estimated PSF");
                        psfSequence.getFirstViewer().getLut().getLutChannel(0).setColorMap(new IceColorMap(),false);
                    }
                });
            }
        };
        workerThread.start();

    }


    private void phaseClicked()
    {

        buildpupil();

        DoubleArray modulus = Double2D.wrap(pupil.getPhi(), new Shape(Nxy, Nxy));
        show(ArrayUtils.roll(modulus),"Phase of the pupil");
    }

    private void modulusClicked()
    {
        /* PSF initialisation */
        buildpupil();
        DoubleArray modulus = Double2D.wrap(pupil.getRho(), new Shape(Nxy, Nxy));
        show(ArrayUtils.roll(modulus),"Modulus of the pupil");
    }

    private void showWeightClicked()
    {
        // Preparing parameters and testing input
        dataSeq = data.getValue();
        dataArray =  sequenceToArray(dataSeq, channel.getValue()).toDouble();
        wgtArray = createWeights(dataArray).toDouble();
        show(wgtArray,"Weight map");
    }


    @Override
    public void stopExecution() {

        run = false;
        if (psfEstimation != null) {
            psfEstimation.stop();
        }
    }


    private void preProcessing(){
        // Preparing parameters and testing input
        dataSeq = data.getValue();
        if (dataSeq == null)
        {
            throwError("An image/sequence of images should be given");
            return;
        }

        dataArray =  sequenceToArray(dataSeq, channel.getValue());
        dataShape = dataArray.getShape();



        if(cursequence==null){
            cursequence = new Sequence("Current Iterate");
            cursequence.copyMetaDataFrom(dataSeq, false);
        }

        Sequence restartSeq = restart.getValue();
        if (restart.getValue() != null && restartSeq != null){
            objArray =  sequenceToArray(restartSeq, channelRestart.getValue());
            if(debug){
                System.out.println("restart seq:" +restartSeq.getName());
            }
        }else{
            if(debug){
                System.out.println("restart seq is null:");
            }
        }

        if (scale.getValue().length !=3){
            throwError("Pixel scale must have 3 elements");
            return;
        }

        if(singlePrecision.getValue()){
            wgtArray = createWeights(dataArray.toFloat()).toFloat();
        }else{
            wgtArray = createWeights(dataArray.toDouble()).toDouble();
        }
    }

    private void deconv() {
        solver = new EdgePreservingDeconvolution();

        solver.setInitialSolution(objArray);

        if(singlePrecision.getValue()){
            solver.setForceSinglePrecision(true);
        }else{
            solver.setForceSinglePrecision(false);
        }

        if(cursequence==null){
            cursequence = new Sequence("Current Iterate");
            cursequence.copyMetaDataFrom(data.getValue(), false);
        }

        solver.setRelativeTolerance(0.0);
        solver.setUseNewCode(false);
        solver.setObjectShape(outputShape);
        solver.setPSF(psfArray);
        solver.setData(dataArray);
        solver.setWeights(wgtArray);
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
        if (debug){
            System.out.println("Launch it:"+nbIterDeconv.getValue());
        }
        run = true;
        OptimTask task = solver.start();

        while (run) {
            task = solver.getTask();
            if (task == OptimTask.ERROR) {
                System.err.format("Error: %s\n", solver.getReason());
                break;
            }
            if (task == OptimTask.NEW_X || task == OptimTask.FINAL_X) {
                if(showIteration.getValue()){
                    show(ArrayUtils.crop(solver.getSolution(),dataShape),cursequence,"Current mu="+solver.getRegularizationLevel() +"it:"+solver.getIterations());
                }
                if (debug){
                    System.out.println("Cost "+solver.getCost() );
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
        objArray = solver.getBestSolution().asShapedArray();
        show(ArrayUtils.crop(objArray,dataShape),cursequence,"Deconvolved "+ dataSeq.getName() + " mu="+solver.getRegularizationLevel());
        solver = null;
    }


    @Override
    public void declareInput(VarList inputMap) {// FIXME Use subclasses for protocols
        initialize();
        super.declareInput(inputMap);
        inputMap.add("dxy", dxy_nm.getVariable());
        inputMap.add("dz", dz_nm.getVariable());
        inputMap.add("NA", na.getVariable());
        inputMap.add("ni", ni.getVariable());
        inputMap.add("lambda", lambda.getVariable());



        inputMap.add("espilon", epsilon.getVariable());

        inputMap.add("padding xy", paddingSizeXY.getVariable());
        inputMap.add("padding z", paddingSizeZ.getVariable());

        inputMap.add("radial",radial.getVariable());
        inputMap.add("nbAlphaCoef", nbAlphaCoef.getVariable());
        inputMap.add("nbBetaCoef", nbBetaCoef.getVariable());
        inputMap.add("defocusMaxIter", maxIterDefocus.getVariable());
        inputMap.add("phaseMaxIter", maxIterPhase.getVariable());
        inputMap.add("modulusMaxIter", maxIterModulus.getVariable());
        inputMap.add("bDecTotalIteration", totalNbOfBlindDecLoop.getVariable());

        inputMap.add("loadFile", loadFile.getVariable());

        deadPixel.setNoSequenceSelection();
    }
    @Override
    public void declareOutput(VarList outputMap) {
        super.declareOutput(outputMap);
        outputMap.add("outputPSF", outputHeadlessPSF.getVariable());

        outputMap.add("pupilShift",pupilShift.getVariable());
        outputMap.add("phaseCoefs",phaseCoefs.getVariable());
        outputMap.add("modulusCoefs",modulusCoefs.getVariable());
    }


    private void parseCmdLine(){
        String[] args = Icy.getCommandLinePluginArgs();

        loadFile.setValue(new File(args[0]));
        loadParamClicked();
        System.out.println("Load Param... "+args[0]);
        for (int i = 1; i < args.length; i++) {
            switch (args[i]) {
                case "-i":
                    if(i+1 >= args.length)
                        break;
                    if( args[i+1].startsWith("-"))
                        break;

                    System.out.println("load image:" + args[i+1]);
                    data.setValue(Loader.loadSequence(args[i+1], 0, false));

                    dataChanged();
                    if(i+3 >= args.length)
                        break;
                    if(args[i+2].equalsIgnoreCase("-c")){
                        channel.setValue(Integer.parseInt(args[i+3]));
                        i=i+3;
                    }else{
                        i++;
                    }


                    break;

                case "-r":
                    if(i+1 >= args.length)
                        break;
                    if( args[i+1].startsWith("-"))
                        break;
                    System.out.println("load restart:" + args[i+1]);
                    restart.setValue(Loader.loadSequence(args[i+1], 0, false));
                    if(i+3 >= args.length){
                        break;}
                    if(args[i+2].equalsIgnoreCase("-c")){
                        System.out.println("channel restart:" + Integer.parseInt(args[i+3]));
                        channelRestart.setValue(Integer.parseInt(args[i+3]));
                        i=i+3;
                    }else{
                        i++;
                    }


                    break;
                case "-o":
                    if(i+1 >= args.length)
                        break;
                    if( args[i+1].startsWith("-"))
                        break;
                    outputPath = args[i+1];
                    i++;
                    break;
                case "-p":
                    if(i+1 >= args.length)
                        break;
                    if( args[i+1].startsWith("-"))
                        break;
                    psfPath = args[i+1];
                    i++;
                    break;
                case "-s":
                    if(i+1 >= args.length)
                        break;
                    if( args[i+1].startsWith("-"))
                        break;
                    saveFile.setValue(new File( args[i+1]));
                    i++;
                    break;
                case "-badpix":
                    if(i+1 >= args.length)
                        break;
                    if( args[i+1].startsWith("-"))
                        break;
                    deadPixel.setValue(Loader.loadSequence(args[i+1], 0, false));
                    i++;
                    break;
                case "-wghtmap":
                    if(i+1 >= args.length)
                        break;
                    if( args[i+1].startsWith("-"))
                        break;
                    weights.setValue(Loader.loadSequence(args[i+1], 0, false));
                    i++;
                    break;

                default:
                    System.out.println("Wrong command line");
                    System.out.println("-i input data file");
                    System.out.println("-r restart file");
                    System.out.println("-o deconvolved output file");
                    System.out.println("-p psf output file");
                    System.out.println("-s parameter output file");
                    System.out.println("-badpix bad pixels file");
                    System.out.println("-wghtmap weight or variance map file");
                    break;
            }
        }

    }
    /**
     *
     */
    @Override
    protected void dataChanged() {
        super.dataChanged();
        pupil=null;
    }
}


