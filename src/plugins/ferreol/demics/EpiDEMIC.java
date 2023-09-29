/*  Copyright (C) 2017  Ferreol Soulez ferreol.soulez@univ-lyon1.fr
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  */

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
import icy.util.OMEUtil;
import microTiPi.epifluorescence.WideFieldModel;
import microTiPi.microUtils.BlindDeconvJob;
import microTiPi.microscopy.PSF_Estimation;
import mitiv.weights.weightsFromModel;
import mitiv.array.ArrayUtils;
import mitiv.array.Double2D;
import mitiv.array.DoubleArray;
import mitiv.base.Shape;
import mitiv.conv.WeightedConvolutionCost;
import mitiv.cost.DifferentiableCostFunction;
import mitiv.cost.HyperbolicTotalVariation;
import mitiv.jobs.DeconvolutionJob;
import mitiv.utils.Histogram;
import ome.units.UNITS;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzButton;
import plugins.adufour.ezplug.EzGroup;
import plugins.adufour.ezplug.EzLabel;
import plugins.adufour.ezplug.EzPanel;
import plugins.adufour.ezplug.EzTabs;
import plugins.adufour.ezplug.EzTabs.TabPlacement;
import plugins.adufour.ezplug.EzVar;
import plugins.adufour.ezplug.EzVarBoolean;
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
public class EpiDEMIC extends DEMICSPlug {


    /***************************************************/
    /**         Plugin interface variables            **/
    /***************************************************/
    private EzTabs tabbedPane;              // The interface is composed of several tabs

    /** data tab: **/
    private EzPanel         dataPanel;      // data tab
    private EzButton        saveMetaData, showPSF;
    private EzVarDouble     dxy_nm;

    /** Noise model tab **/
    private EzPanel         noisePanel; // Panel with noise parameters



    /** deconvolution tab: **/
    private EzPanel         deconvPanel;
    private EzVarDouble     epsilon; // deconvolution hyper parameters; mu = 10^(logmu)

    protected  int          Nxy=32; // Output (padded sequence size)


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
        super.initialize();

        tabbedPane = new EzTabs("BlindTabs", TabPlacement.TOP);

        /****************************************************/
        /**                    IMAGE TAB                   **/
        /****************************************************/

        dataPanel = new EzPanel("Step 1: Data");
        saveMetaData = new EzButton("Save metadata", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if (debug) {
                    System.out.println("Saving metadata");
                }
                updateMetaData();
            }
        });


        paddingSizeXY = new EzVarInteger("padding in x and y:",0, Integer.MAX_VALUE,1);

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



        dataEV.addVarChangeListener(new EzVarListener<Sequence>() {
            @Override
            public void variableChanged(EzVar<Sequence> source,
                    Sequence newValue) {
                if (debug) {
                    System.out.println("Seq ch..."+newValue);
                }
                dataChanged() ;
                if (sizeZ == 1) {
                    startDecButton.setEnabled(false);
                    startBlind.setEnabled(false);
                    throwError("Input channel must be 3D");
                    return;
                }
                if(dataSeq!=null){
                    startBlind.setEnabled(true);
                    if (dataSeq.getPixelSizeX()!=dataSeq.getPixelSizeY()) {
                        startDecButton.setEnabled(false);
                        startBlind.setEnabled(false);
                        throwError("Input channel must have the same scale in X and Y");
                        return;
                    }
                    dxy_nm.setValue( dataSeq.getPixelSizeX()*1E3);
                    dz_nm.setValue( dataSeq.getPixelSizeZ()*1E3);
                    scale.setValue(new double[]{1.0 ,1.0, dxy_nm.getValue()/ dz_nm.getValue() } );

                    try {
                        lambda.setValue( metDat.getChannelEmissionWavelength(0,channelEV.getValue()).value(UNITS.NANOMETER).doubleValue());
                    } catch(Exception e){
                        System.out.println("Failed to get some wavelength metadatas, will use default values ");
                        lambda.setValue(500.0);
                    }
                    try {
                        na.setValue( metDat.getObjectiveLensNA(0, 0));
                    } catch(Exception e){
                        System.out.println("Failed to get na metadatas, will use default values ");
                        na.setValue(1.4);
                    }

                    try {
                        if (metDat.getObjectiveSettingsRefractiveIndex(0)!=null)
                            ni.setValue(metDat.getObjectiveSettingsRefractiveIndex(0) );
                        else {
                            System.out.println("Failed to get refractive index from metadata, will use default values ");
                            ni.setValue(1.518);
                        }
                    } catch(Exception e){
                        System.out.println("Failed to get refractive index from metadata, will use default values ");
                        ni.setValue(1.518);
                    }
                    if (debug) {
                        System.out.println("Seq changed:" + sizeX + "  "+ Nxy);
                    }
                }else{
                    startDecButton.setEnabled(false);
                    startBlind.setEnabled(false);
                }

            }

        });


        EzVarListener<Double> metaActionListener = new EzVarListener<Double>() {
            @Override
            public void variableChanged(EzVar<Double> source, Double newValue) {
                scale.setValue(new double[]{1.0 ,1.0,  dz_nm.getValue()/dxy_nm.getValue() } );
                pupilShift.setValue(new double[] { 0., 0.});
                pupil=null;
            };
        };



        dxy_nm = new EzVarDouble("dxy(nm):",64.5,0., Double.MAX_VALUE,1.);
        dxy_nm.addVarChangeListener(metaActionListener);
        dx_nm = new EzVarDouble("dx(nm):",64.5,0., Double.MAX_VALUE,1.);
        dx_nm.addVarChangeListener(metaActionListener);
        dy_nm = new EzVarDouble("dy(nm):",64.5,0., Double.MAX_VALUE,1.);
        dy_nm.addVarChangeListener(metaActionListener);
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



        ezPaddingGroup = new EzGroup("Padding",paddingSizeXY, paddingSizeZ);
        ezPaddingGroup.setFoldedState(true);



        /****************************************************/
        /**                    NOISE  TAB                  **/
        /****************************************************/
        noisePanel = new EzPanel("Step 2: Noise model");

        /****************************************************/
        /**                    DECONV TAB                  **/
        /****************************************************/
        deconvPanel = new EzPanel("Step 3: Deconvolution");
        epsilon = new EzVarDouble("Threshold level:",1E-2,0.,Double.MAX_VALUE,1.0);
        singlePrecision.addVarChangeListener(new EzVarListener<Boolean>() {
            @Override
            public void variableChanged(EzVar<Boolean> source, Boolean newValue) {
                pupil=null;
                psfEstimation=null;
                cursequence =null;

            }
        });
        docDec = new EzLabel("Launch a deconvolution using the current PSF", Color.red);
        startDecButton = new EzButton("Start Deconvolution", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                launchClicked(false);
            }
        });


        ezDeconvolutionGroup = new EzGroup("Expert  parameters",epsilon,scale,positivityEV,singlePrecision);//, showFullObjectButton);
        ezDeconvolutionGroup.setFoldedState(true);

        /****************************************************/
        /**                      BDEC TAB                  **/
        /****************************************************/
        bdecPanel = new EzPanel("Step 4: Blind deconvolution");
        pupilShift = new EzVarDoubleArrayNative("pupilShift", new double[][] { { 0.0, 0.0} }, false);
        pupilShift.setVisible(false);
        phaseCoefs = new EzVarDoubleArrayNative("phase coefs",null , false);
        phaseCoefs.setVisible(false);
        modulusCoefs = new EzVarDoubleArrayNative("modulusCoefs", new double[][] { {1.0 } },0, false);
        modulusCoefs.setVisible(false);

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
                double[] tmp = new double[Math.max(1,Integer.parseInt(nbBetaCoef.getValue()))];
                tmp[0] = 1;
                modulusCoefs.setValue(tmp);
                pupilShift.setValue(new double[] { 0., 0.});
                if (metDat!=null)
                    ni.setValue(1.518); //FIXME
                else
                    ni.setValue(1.518);
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
                fSeq.copyMetaDataFrom(dataEV.getValue(), false);
                fSeq.setMetaData(OMEUtil.createOMEXMLMetadata(dataEV.getValue().getOMEXMLMetadata()));
                if(objArray != null){
                    IcyImager.show(objArray,fSeq,"Deconvolved "+ dataEV.getValue().getName() + " with padding. mu " + mu.getValue(),isHeadLess());
                }else {
                    IcyImager.show(ArrayUtils.extract(dataArray, outputShape),fSeq,"Deconvolved "+ dataEV.getValue().getName() + "with padding. mu="+ mu.getValue(),isHeadLess() );
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
        dataEV.setToolTipText(ToolTipText.sequenceImage);
        weightsSeq.setToolTipText(ToolTipText.sequenceWeigth);
        weightsMethod.setToolTipText(ToolTipText.sequenceWeigth);
        badpixMap.setToolTipText(ToolTipText.sequencePixel);

        channelEV.setToolTipText(ToolTipText.textCanal);

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

        restartEV.setToolTipText(ToolTipText.booleanRestart);
        positivityEV.setToolTipText(ToolTipText.booleanPositivity);
        showFullObjectButton.setToolTipText(ToolTipText.booleanCrop);
        showFullObject2.setToolTipText(ToolTipText.booleanCrop);
        debuggingPanel.setToolTipText(ToolTipText.textOutput);



        /**** IMAGE ****/
        dataPanel.add(dataEV);
        dataPanel.add(channelEV);
        dataPanel.add(dataSizeTxt);
        dataPanel.add(ezPaddingGroup);
        dataPanel.add(outputSizeTxt);
        dataPanel.add(dxy_nm);
        dataPanel.add(dz_nm);


        dataPanel.add(na);
        dataPanel.add(ni);
        dataPanel.add(lambda);

        //    dataPanel.add(ezWeightingGroup);

        dataPanel.add(saveMetaData);
        dataPanel.add(showPSF);
        dataPanel.add(loadFile);
        dataPanel.add(loadParam);
        tabbedPane.add(dataPanel);

        /**** Noise model ****/
        noisePanel.add(new EzLabel("The noise is supposed Gaussian with variance given by:"
                + "\n"));
        noisePanel.add(weightsMethod);
        noisePanel.add(weightsSeq);
        noisePanel.add(gain);
        noisePanel.add(noise);
        noisePanel.add(badpixMap);
        noisePanel.add(showWeightButton);
        tabbedPane.add(noisePanel);


        /**** Deconv ****/

        deconvPanel.add(logmu);
        deconvPanel.add(mu);
        deconvPanel.add(nbIterDeconv);
        deconvPanel.add(restartEV);
        deconvPanel.add(channelRestartEV);
        deconvPanel.add(ezDeconvolutionGroup);

        deconvPanel.add(docDec);
        deconvPanel.add(startDecButton);
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

        if (!isHeadLess()) {
            outputSizeTxt.setEnabled(false);
            dataSizeTxt.setEnabled(false);
            mu.setEnabled(false);
        }else{
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
                            restartEV.setValue(cursequence);
                            channelRestartEV.setValue(0);
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
    @Override
    protected void loadParamClicked() {
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

    @Override
    protected void saveParamClicked() {
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
            dataChanged();
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
            startDecButton.setText("Emergency stop");
            buildpupil();
            if (debug|| isHeadLess()) {
                System.out.println("-------------IMAGE-------------------");
                System.out.println("File: "+dataEV.getValue());
                System.out.println("Canal: "+channelEV.getValue());
                System.out.println("image size: "+ dataSizeTxt.getValue());
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
                System.out.println("Weights: "+weightsSeq.getValue());
                System.out.println("Gain: "+gain.getValue());
                System.out.println("Noise: "+noise.getValue());
                System.out.println("deadPix: "+badpixMap.getValue());
                System.out.println("--------------DECONV------------------");

                System.out.println("zeroPad xy: "+paddingSizeXY.getValue());
                System.out.println("zeroPad z: "+paddingSizeZ.getValue());
                System.out.println("nbIter: "+nbIterDeconv.getValue());
                System.out.println("Restart: "+restartEV.getValue());
                System.out.println("Positivity: "+positivityEV.getValue());
                System.out.println("--------------BDEC------------------");
                System.out.println("output size: "+ outputSizeTxt.getValue());
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


                psfEstimation = new PSF_Estimation(pupil);

                psfArray = ArrayUtils.roll( pupil.getPsf() );
                preProcessing();

                psfEstimation.setData(ArrayUtils.pad(dataArray,outputShape));
                psfEstimation.setWeight(  ArrayUtils.pad(wgtArray,outputShape));

                psfEstimation.enablePositivity(false);
                psfEstimation.setAbsoluteTolerance(0.0);

                int[] bMaxIter = {maxIterDefocus.getValue(),maxIterPhase.getValue(), maxIterModulus.getValue()};


                bdec = new BlindDeconvJob(totalNbOfBlindDecLoop.getValue(), pupil.getParametersFlags(), bMaxIter, psfEstimation ,deconvolver, wghtUpdt,debug );

                objArray = bdec.blindDeconv(objArray);
                if(wghtUpdt!=null){
                    wgtArray = wghtUpdt.getWeights();
                    gain.setValue(((weightsFromModel) wghtUpdt).getAlpha());
                    noise.setValue(Math.sqrt(((weightsFromModel) wghtUpdt).getBeta())/gain.getValue());
                }
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
                    restartEV.setValue(cursequence);
                    channelRestartEV.setValue(0);
                    if (isHeadLess()) {
                        if(outputHeadlessImage==null){
                            outputHeadlessImage = new EzVarSequence("Output Image");
                        }
                        if(outputHeadlessPSF==null){
                            outputHeadlessPSF = new EzVarSequence("Output PSF");
                        }

                        if (outputHeadlessWght==null) {
                            outputHeadlessWght = new EzVarSequence("Computed weightsSeq");
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
            startDecButton.setText("Start Deconvolution");
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

            dataEV.setEnabled(flag);
            channelEV.setEnabled(flag);
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
        dataSeq = dataEV.getValue();
        if (singlePrecision.getValue()) {
            dataArray =  sequenceToArray(dataSeq, channelEV.getValue()).toFloat();
        }else {
            dataArray =  sequenceToArray(dataSeq, channelEV.getValue()).toDouble();
        }
        createWeights(true);
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
        dataSeq = dataEV.getValue();
        if (dataSeq == null)
        {
            throwError("No image/sequence");
            return;
        }


        if (singlePrecision.getValue()) {
            dataArray =  sequenceToArray(dataSeq, channelEV.getValue()).toFloat();
        }else {
            dataArray =  sequenceToArray(dataSeq, channelEV.getValue()).toDouble();
        }
        dataShape = dataArray.getShape();


        Sequence restartSeq = restartEV.getValue();
        if (  restartSeq != null){
            objArray =  sequenceToArray(restartSeq, channelRestartEV.getValue());
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

        createWeights(true);
        if (weightsMethod.getValue() == weightOptions[4]) {
            wghtUpdt = new weightsFromModel( dataArray, badpixArray);
        }

        if(cursequence==null){
            cursequence = new Sequence("Current Iterate");
            cursequence.copyMetaDataFrom(dataEV.getValue(), false);
        }
        IcyImager curImager = new IcyImager(cursequence, isHeadLess());

        DeconvHook dHook = new DeconvHook(curImager, dataShape,null, debug);
        DeconvHook dHookfinal = new DeconvHook(curImager, dataShape,"Deconvolved "+dataSeq.getName(), debug);

        buildVectorSpaces();

        DifferentiableCostFunction fprior = new HyperbolicTotalVariation(objectSpace, epsilon.getValue(), scale.getValue());
        WeightedConvolutionCost fdata =  WeightedConvolutionCost.build( objectSpace, dataSpace);
        fdata.setData(dataArray);
        fdata.setWeights(wgtArray,true);
        fdata.setPSF(psfArray);
        deconvolver  = new DeconvolutionJob( fdata,  mu.getValue(),fprior,  positivityEV.getValue(),nbIterDeconv.getValue(),  dHook,  dHookfinal);
        objArray = ArrayUtils.extract(objArray, outputShape, 0.); //Padding to the right size


    }

    private void deconv() {
        if (debug){
            System.out.println("Launch it:"+nbIterDeconv.getValue());
        }
        objArray = deconvolver.deconv(objArray);
        if(wghtUpdt!=null) {
            wghtUpdt.update(deconvolver);
            wgtArray = wghtUpdt.getWeights();
            gain.setValue(((weightsFromModel) wghtUpdt).getAlpha());
            noise.setValue(Math.sqrt(((weightsFromModel) wghtUpdt).getBeta())/gain.getValue());
        }
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

        badpixMap.setNoSequenceSelection();
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
                    dataEV.setValue(Loader.loadSequence(args[i+1], 0, false));

                    dataChanged();
                    if(i+3 >= args.length)
                        break;
                    if(args[i+2].equalsIgnoreCase("-c")){
                        channelEV.setValue(Integer.parseInt(args[i+3]));
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
                    restartEV.setValue(Loader.loadSequence(args[i+1], 0, false));
                    if(i+3 >= args.length){
                        break;}
                    if(args[i+2].equalsIgnoreCase("-c")){
                        System.out.println("channel restart:" + Integer.parseInt(args[i+3]));
                        channelRestartEV.setValue(Integer.parseInt(args[i+3]));
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
                    badpixMap.setValue(Loader.loadSequence(args[i+1], 0, false));
                    i++;
                    break;
                case "-wghtmap":
                    if(i+1 >= args.length)
                        break;
                    if( args[i+1].startsWith("-"))
                        break;
                    weightsSeq.setValue(Loader.loadSequence(args[i+1], 0, false));
                    i++;
                    break;

                default:
                    System.out.println("Wrong command line");
                    System.out.println("-i input channel file");
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
        if(dataSeq!=null){
            epsilon.setValue((new Histogram(sequenceToArray(dataSeq, channelEV.getValue()))).getMaximumValue()/1000 );
        };
        badpixArray = null;
        psfEstimation=null;
        cursequence =null;
        pupil=null;
        deconvolver = null;
    }

}
