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
import icy.plugin.PluginDescriptor;
import icy.plugin.PluginInstaller;
import icy.plugin.PluginRepositoryLoader;
import icy.plugin.PluginUpdater;
import icy.sequence.Sequence;
import icy.system.thread.ThreadUtil;
import loci.formats.ome.OMEXMLMetadata;
import microTiPi.epifluorescence.WideFieldModel;
import microTiPi.microUtils.BlindDeconvJob;
import microTiPi.microscopy.MicroscopeMetadata;
import microTiPi.microscopy.PSF_Estimation;
import mitiv.array.ArrayUtils;
import mitiv.array.Double2D;
import mitiv.array.DoubleArray;
import mitiv.array.ShapedArray;
import mitiv.base.Shape;
import mitiv.base.mapping.DoubleFunction;
import mitiv.jobs.DeconvolutionJob;
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
import plugins.mitiv.io.DeconvHook;
import plugins.mitiv.io.IcyImager;

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

    protected MicroscopeMetadata meta = null; // metadata of the data

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
    private EzButton        startBlind,  showFullObject2, resetPSF;
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
    private DeconvolutionJob deconvolver ;


    // Main arrays for the psf estimation
    //  private boolean runBdec;
    private WideFieldModel pupil=null;


    /*********************************/
    /**            DEBUG            **/
    /*********************************/
    private boolean debug = false;      // Show psf steps
    private EzPanel  debuggingPanel;
    private EzVarText resultCostPrior, resultDefocus, resultPhase, resultModulus;


    private String psfPath=null;

    private BlindDeconvJob bdec;

    private ShapedArray badArray=null;




    @Override
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
        Nx  = Ny = Nxy;
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

        // REMOVE old MitivDeconvolution plugin
        PluginRepositoryLoader.waitLoaded();
        // get  plugin descriptor
        PluginDescriptor desc = PluginRepositoryLoader.getPlugin("plugins.mitiv.deconv.MitivDeconvolution");
        // install  plugin
        if(desc != null){
            if( desc.isInstalled()){
                if (!isHeadLess()) {
                    new AnnounceFrame("Removing the now useless MitivDeconvolution plugin.");
                }
                PluginInstaller.desinstall(desc, false, false);
                while (PluginUpdater.isCheckingForUpdate() ||  PluginInstaller.isProcessing() || PluginInstaller.isInstalling() || PluginInstaller.isDesinstalling())
                    ThreadUtil.sleep(1);
            }
        }
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
                    startDec.setEnabled(true);
                    startBlind.setEnabled(true);
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
                }else{
                    startDec.setEnabled(false);
                    startBlind.setEnabled(false);
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
        deadPixel.addVarChangeListener(new EzVarListener<Sequence>() {

            @Override
            public void variableChanged(EzVar<Sequence> source, Sequence newValue) {
                if (newValue!=null){
                    badArray = sequenceToArray(newValue);
                }else{
                    badArray = null;
                }
            }
        });
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
                launchClicked(false);
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
                    IcyImager.show(objArray,fSeq,"Deconvolved "+ data.getValue().getName() + " with padding. mu " + mu.getValue(),isHeadLess() );
                }else {
                    IcyImager.show(ArrayUtils.extract(dataArray, outputShape),fSeq,"Deconvolved "+ data.getValue().getName() + "with padding. mu="+ mu.getValue(),isHeadLess() );
                }
            }
        });


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
                phaseCoefs.setValue(new double[Integer.parseInt(nbAlphaCoef.getValue())]);
                double[] tmp = new double[Integer.parseInt(nbBetaCoef.getValue())];
                tmp[0] = 1;
                modulusCoefs.setValue(tmp);
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
                    nbBetaCoef.setDefaultValues(  nBetaOptions,0, false );
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
                launchClicked(true);

            }
        });



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
                    IcyImager.show(objArray,fSeq,"Deconvolved "+ data.getValue().getName() + " with padding. mu " + mu.getValue(),isHeadLess());
                }else {
                    IcyImager.show(ArrayUtils.extract(dataArray, outputShape),fSeq,"Deconvolved "+ data.getValue().getName() + "with padding. mu="+ mu.getValue(),isHeadLess() );
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
        tabbedPane.add(deconvPanel);

        /**** BDec ****/
        bdecPanel.add(resetPSF);
        bdecPanel.add(totalNbOfBlindDecLoop);
        bdecPanel.add(ezBlindDeconvolutionGroup);
        bdecPanel.add(docBlind);
        bdecPanel.add(startBlind);
        bdecPanel.add(visuPSF);
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
     * @param flag
     */
    private void launchClicked(final boolean flag) {

        if ( deconvolver!=null && deconvolver.isRunning()){
            stopExecution();
        } else if ( bdec!=null && bdec.isRunning()){
            stopExecution();
        }else{
            Thread  workerThread = new Thread() {
                @Override
                public void run() {
                    launch(flag);

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


    }





    /**
     * Load the parameter file and perform parameter update
     */
    private void loadParamClicked() {
        File loadName = loadFile.getValue();

        if ( !isHeadLess()){
            new AnnounceFrame("Loading deconvolution parameters from "+ loadName.getAbsolutePath().toString(),3);
        }
        this.loadParameters(loadFile.getValue());
        loadFile.setValue(loadName);    // FIX saving file name erasing during load

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
            if ( !isHeadLess()){
                new AnnounceFrame("Saving deconvolution parameters in "+(pathName.getAbsolutePath().toString()),3);
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
        launch(true);

        long stopTime = System.currentTimeMillis();
        long elapsedTime = stopTime - startTime;
        System.out.println("time: "+elapsedTime);
    }

    private void launch(boolean runBdec) {
        try {
            startBlind.setText("Emergency stop");
            startDec.setText("Emergency stop");
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


            /*---------------------------------------*/
            /*            OPTIMISATION               */
            /*---------------------------------------*/

            if (runBdec) {

                int nbAlpha = Integer.parseInt(nbAlphaCoef.getValue());
                int nbBeta = Integer.parseInt(nbBetaCoef.getValue());

                if (nbAlpha==0){
                    maxIterPhase.setValue(0);
                }else{
                    pupil.setNPhase(nbAlpha);
                }

                if (nbBeta==0){
                    maxIterModulus.setValue(0);
                }else{
                    pupil.setNModulus(nbBeta);
                }

                preProcessing();

                psfEstimation = new PSF_Estimation(pupil);

                psfEstimation.setWeight(  ArrayUtils.pad(wgtArray,outputShape));
                psfEstimation.setData(ArrayUtils.pad(dataArray,outputShape));

                psfEstimation.enablePositivity(false);
                psfEstimation.setAbsoluteTolerance(0.0);

                int[] bMaxIter = {maxIterDefocus.getValue(),maxIterPhase.getValue(), maxIterModulus.getValue()};


                bdec = new BlindDeconvJob(totalNbOfBlindDecLoop.getValue(), pupil.getParametersFlags(), bMaxIter, psfEstimation ,deconvolver, debug );


                objArray = bdec.blindDeconv(objArray);

                if(maxIterDefocus.getValue()>0){
                    ni.setValue(((WideFieldModel) psfEstimation.getModel()).getNi());
                    pupilShift.setValue(((WideFieldModel) psfEstimation.getModel()).getPupilShift());
                }
                if(maxIterPhase.getValue()>0){
                    phaseCoefs.setValue(((WideFieldModel) psfEstimation.getModel()).getPhaseCoefs().getData());
                }
                if(maxIterModulus.getValue()>0){
                    modulusCoefs.setValue(((WideFieldModel) psfEstimation.getModel()).getModulusCoefs().getData());
                }
                pupil =((WideFieldModel) psfEstimation.getModel());

                bdec = null;
            } else {
                psfArray = ArrayUtils.roll( pupil.getPsf() );
                pupil.freeMem();
                preProcessing();
                deconv();
            }
            if(pupil!=null)
                pupil.freeMem();// TODO free more memory

            SwingUtilities.invokeLater(new Runnable() {
                @Override
                public void run() {
                    enableVars(true);
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
                        psfSequence =   arrayToSequence( ArrayUtils.roll(pupil.getPsf()), psfSequence);


                        outputHeadlessPSF.setValue(psfSequence);
                        outputHeadlessImage.setValue(cursequence);
                        outputHeadlessWght.setValue(arrayToSequence(wgtArray));

                        if(outputPath!=null){
                            IcyImager.save(cursequence, outputPath);
                        }

                        if(psfPath!=null){
                            IcyImager.save(psfSequence, psfPath);
                        }

                        if(saveFile.getValue()!=null){
                            saveParamClicked();
                        }
                    }
                }
            });
        } catch (IllegalArgumentException e) {
            if(!isHeadLess()){
                new AnnounceFrame("Oops, Error: "+ e.getMessage());
            }
            enableVars(true);

            if (debug) {
                e.printStackTrace();
            }
        } finally {
            startBlind.setText("Guess PSF");
            startDec.setText("Start Deconvolution");
        }
    }


    /**
     * If false disable some variables in the interface
     * @param flag
     *
     */
    private void enableVars(boolean flag) {

        if (!isHeadLess()) {
            nbAlphaCoef.setEnabled(flag);
            nbBetaCoef.setEnabled(flag);
            radial.setEnabled(flag);

            na.setEnabled(flag);
            lambda.setEnabled(flag);
            dxy_nm.setEnabled(flag);
            dz_nm.setEnabled(flag);
            ni.setEnabled(flag);

            data.setEnabled(flag);
            channel.setEnabled(flag);
            paddingSizeXY.setEnabled(flag);
            paddingSizeZ.setEnabled(flag);

            singlePrecision.setEnabled(flag);
            loadParam.setEnabled(flag);
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
                if (pupil==null){
                    buildpupil();
                }

                SwingUtilities.invokeLater(new Runnable() {
                    @Override
                    public void run() {
                        Sequence psfSequence = new Sequence();
                        if(dataSeq!=null){
                            psfSequence.copyMetaDataFrom(dataSeq, false);
                        }

                        if ( bdec!=null && bdec.isRunning()){
                            IcyImager.show( bdec.getPsf(),psfSequence,"Estimated PSF",isHeadLess());
                        }else{
                            IcyImager.show(ArrayUtils.roll(pupil.getPsf()),psfSequence,"Estimated PSF",isHeadLess());
                        }
                        psfSequence.getFirstViewer().getLut().getLutChannel(0).setColorMap(new IceColorMap(),false);
                    }
                });
            }
        };
        workerThread.start();

    }


    private void phaseClicked()
    {
        DoubleArray  phase;
        if (pupil==null){
            buildpupil();
        }
        if ( bdec!=null && bdec.isRunning()){
            phase = Double2D.wrap(((WideFieldModel) bdec.getPupil()).getPhi(), new Shape(Nxy, Nxy));
        }else{
            phase = Double2D.wrap(pupil.getPhi(), new Shape(Nxy, Nxy));
        }
        IcyImager.show(ArrayUtils.roll(phase),null,"Phase of the pupil",false);
    }

    private void modulusClicked()
    {
        DoubleArray modulus;
        if (pupil==null){
            buildpupil();
        }

        if ( bdec!=null && bdec.isRunning()){
            modulus = Double2D.wrap(((WideFieldModel) bdec.getPupil()).getRho(), new Shape(Nxy, Nxy));
        }else{
            modulus = Double2D.wrap(pupil.getRho(), new Shape(Nxy, Nxy));
        }
        IcyImager.show(ArrayUtils.roll(modulus),null,"Modulus of the pupil",false);
    }

    private void showWeightClicked()
    {
        // Preparing parameters and testing input
        dataSeq = data.getValue();
        dataArray =  sequenceToArray(dataSeq, channel.getValue()).toDouble();
        wgtArray = createWeights(dataArray,badArray).toDouble();
        IcyImager.show(wgtArray,null,"Weight map",false);
    }


    @Override
    public void stopExecution() {
        if (deconvolver!=null)
            deconvolver.abort();
        //  run = false;
        if(bdec!=null)
            bdec.abort();
    }


    private void preProcessing(){
        // Preparing parameters and testing input
        dataSeq = data.getValue();
        if (dataSeq == null)
        {
            throwError("No image/sequence");
            return;
        }

        dataArray =  sequenceToArray(dataSeq, channel.getValue());
        dataShape = dataArray.getShape();
        if (deadPixel.getValue() ==null){
            badArray = null;
            if (dataSeq.getChannelMax( channel.getValue())>=( dataSeq.getChannelTypeMax(channel.getValue())-1)){
                class SaturationFunc implements DoubleFunction{
                    double sat;
                    public  SaturationFunc(double sat){
                        this.sat = sat;
                    }
                    @Override
                    public double apply(double arg) {
                        if  (arg>= this.sat){
                            return 1.;
                        }else{
                            return 0.;
                        }
                    }

                }

                badArray = dataArray.copy().toDouble();
                ((DoubleArray) badArray).map(new SaturationFunc(dataSeq.getChannelTypeMax(channel.getValue())));
                badArray = badArray.toByte();
                if (!isHeadLess()){
                    new AnnounceFrame("Warning, saturated pixel detected, accounting them as dead pixels", "show", new Runnable() {
                        @Override
                        public void run() {
                            Sequence deadSequence = new Sequence("Saturations map");
                            deadSequence.copyMetaDataFrom(dataSeq, false);
                            IcyImager.show(badArray, deadSequence, "saturations map", isHeadLess());
                        }
                    }, 10);
                }
            }
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
            wgtArray = createWeights(dataArray.toFloat(),badArray).toFloat();
        }else{
            wgtArray = createWeights(dataArray.toDouble(),badArray).toDouble();
        }

        if (scale.getValue().length !=3){
            throwError("Pixel scale must have 3 elements");
            return;
        }

        if(cursequence==null){
            cursequence = new Sequence("Current Iterate");
            cursequence.copyMetaDataFrom(data.getValue(), false);
        }
        IcyImager curImager = new IcyImager(cursequence, isHeadLess());

        DeconvHook dHook = new DeconvHook(curImager, dataShape,null, debug);
        DeconvHook dHookfinal = new DeconvHook(curImager, dataShape,"Deconvolved "+dataSeq.getName(), debug);
        deconvolver = new DeconvolutionJob(dataArray, psfArray, wgtArray, outputShape, mu.getValue(), epsilon.getValue(), scale.getValue(), positivity.getValue(), singlePrecision.getValue(), nbIterDeconv.getValue(), dHook , dHookfinal);

    }

    private void deconv() {
        if (debug){
            System.out.println("Launch it:"+nbIterDeconv.getValue());
        }
        objArray = deconvolver.deconv(objArray);

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

        badArray = null;
        psfEstimation=null;
        cursequence =null;
        pupil=null;
        deconvolver = null;
    }

    /**
     * Here we get the informations given by the users but not all.
     * In fact we trust only a few data that we know that are given by Icy.
     * Else we are trying to keep them for the next run.
     *
     * Remember: if users may lie, they will !
     *
     * @param seq
     * @return
     */
    protected MicroscopeMetadata getMetaData(Sequence seq){ //FIXME Should be elsewhere
        OMEXMLMetadata metDat = seq.getMetadata();
        if (meta == null) {
            meta = new MicroscopeMetadata();
            if (metDat.getInstrumentCount() > 0) {
                try {
                    meta.na      = metDat.getObjectiveLensNA(0, 0);
                    //meta.lambda  = metDat.getChannelEmissionWavelength(0, 0).getValue().doubleValue()*1E6;  //I suppose the value I will get is in um
                } catch(Exception e){
                    System.out.println("Failed to get some metadatas, will use default values for na, lambda");
                }
            }
        }
        //If no instrument found, at least we have the right image size
        meta.nxy     = seq.getSizeX(); //We suppose X and Y equal
        meta.nz      = seq.getSizeZ();
        meta.dxy     = seq.getPixelSizeX()*1E3;
        meta.dz      = seq.getPixelSizeZ()*1E3;
        meta.na      = na.getValue();
        meta.lambda  = lambda.getValue();
        meta.ni      = ni.getValue();
        return meta;
    }

}


