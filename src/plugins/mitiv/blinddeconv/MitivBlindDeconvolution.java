package plugins.mitiv.blinddeconv;

import icy.gui.frame.progress.AnnounceFrame;
import icy.image.IcyBufferedImage;
import icy.plugin.interface_.PluginBundled;
import icy.sequence.MetaDataUtil;
import icy.sequence.Sequence;
import icy.type.collection.array.DynamicArray.Int;
import icy.util.OMEUtil;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Arrays;

import javax.swing.SwingUtilities;

import loci.common.services.ServiceException;
import loci.formats.ome.OMEXMLMetadata;
import loci.formats.ome.OMEXMLMetadataImpl;
import mitiv.array.ArrayFactory;
import mitiv.array.ArrayUtils;
import mitiv.array.Double2D;
import mitiv.array.Double3D;
import mitiv.array.DoubleArray;
import mitiv.array.FloatArray;
import mitiv.array.ShapedArray;
import mitiv.base.Shape;
import mitiv.base.Traits;
import mitiv.invpb.ReconstructionJob;
import mitiv.invpb.ReconstructionViewer;
import mitiv.linalg.shaped.DoubleShapedVector;
import mitiv.linalg.shaped.DoubleShapedVectorSpace;
import mitiv.microscopy.WideFieldModel;
import mitiv.microscopy.PSF_Estimation;
import mitiv.utils.FFTUtils;
import mitiv.utils.MathUtils;
import mitiv.utils.WeightFactory;
import mitiv.utils.reconstruction.ReconstructionThread;
import mitiv.utils.reconstruction.ReconstructionThreadToken;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzButton;
import plugins.adufour.ezplug.EzGroup;
import plugins.adufour.ezplug.EzLabel;
import plugins.adufour.ezplug.EzPanel;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;
import plugins.adufour.ezplug.EzTabs;
import plugins.adufour.ezplug.EzVar;
import plugins.adufour.ezplug.EzTabs.TabPlacement;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarChannel;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarListener;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.ezplug.EzVarText;
import plugins.mitiv.deconv.MitivDeconvolution;
import plugins.mitiv.io.IcyBufferedImageUtils;
import plugins.mitiv.myEzPlug.MyMetadata;
import plugins.mitiv.reconstruction.TotalVariationJobForIcy;

/**
 * MiTivGlobalDeconv is a blind deconvolution tool built on the same basis than
 * MiTivTotalVariation. The blind deconvolution process is trying to guess the PSF 
 * and then the it is a standard deconvolution. 
 * 
 * @author light
 *
 */
public class MitivBlindDeconvolution extends EzPlug implements EzStoppable, Block, PluginBundled {

    /***************************************************/
    /**               Viewer Update result            **/
    /***************************************************/
    public class TvViewer implements ReconstructionViewer{
        @Override
        public void display(ReconstructionJob job) {
            setResult(job);
        }
    }

    /***************************************************/
    /**                  All variables                **/
    /***************************************************/
    private EzVarDouble dxy_nm, dz_nm, na, lambda, ni;
    private EzVarInteger  nbIteration, bDecTotalIteration,defocusMaxIter,phaseMaxIter,modulusMaxIter, paddingSizeXY, paddingSizeZ;
    private WideFieldModel pupil=null;
    // private boolean psfInitFlag = false;
    private EzVarDouble logmu, mu, epsilon, gain, noise;         
    private EzVarSequence image, restart, weights, deadPixel, outputHeadlessImage, outputHeadlessPSF;
    private EzVarText weightsMethod,  nbAlphaCoef, nbBetaCoef;
    private EzVarChannel channel;
    private EzVarBoolean deadPixGiven, positivity,radial; // crop
    private final String[] weightOptions = new String[]{"None","Inverse covariance map","Variance map","Computed variance"}; 
    private final String[] nAlphaOptions = new String[]{"0","1","3","7","8","12","18","19","25","33","34","42","52","53","63","75","76","88","102","103"}; 
    private final String[] nBetaOptions = new String[]{"0","3","4","6","10","11","15","21","22","28","36","37","45","55","56","66","78","79","91","105","106"};
    private final String[] nAlphaOptionsR = new String[]{"0","1","2","3","4","5","6","7","8","9"}; 
    private final String[] nBetaOptionsR = new String[]{"0","1","2","3","4","5","6","7","8","9"};
    private MyMetadata meta = null;     //The image metadata that we will move from one image to another
    private EzButton saveMetaData, showPSF, psfShow2, showWeight, showModulus, showPhase;
    private EzButton startDec, startBlind, stopDec, stopBlind, cropResult, cropResultInDeconv, resetPSF;
    private EzVarText imageSize, outputSize, resultCostPrior, resultDefocus, resultPhase, resultModulus;
    private EzLabel docDec, docBlind;
	private EzVarBoolean expertMode;

    private EzPanel  imageGlob, varianceGlob, deconvGlob, bdecGlob, resultGlob; 
    private EzTabs tabbedPane;

    private double grtol = 0.0;
    private int nbAlpha=0, nbBeta=1;
    private int sizeX=512, sizeY=512, sizeZ=128; // Input sequence sizes 
    private Shape imageShape;
    private  int Nxy=512, Nz=128;			 // Output (padded sequence size)
    private Shape outputShape;
    


    private boolean guessPhase;


    DoubleShapedVectorSpace defocuSpace = null, alphaSpace=null, betaSpace=null;
    DoubleShapedVector defocusVector = null, alphaVector = null, betaVector =null;


    boolean run = true;
    boolean runBdec;

    /*********************************/
    /**            Job              **/
    /*********************************/
    private ReconstructionThreadToken token;
    ReconstructionThread thread;

    /*********************************/
    /**            DEBUG            **/
    /*********************************/
    private boolean debug = false;      // Show psf steps 
    private boolean verbose = false;    // Show some values, need debug to true

    // Global variables for the algorithms
    TotalVariationJobForIcy tvDec;
    PSF_Estimation PSFEstimation;
    //Global variable for the deconvolution
    Sequence sequence; //The reference to the sequence we use to plot 
    private boolean guessModulus;
    Sequence lastSequence; // Just a reference to know the last result
	private EzGroup show;


    /*********************************/
    /**       Utils functions       **/
    /*********************************/

    private static void throwError(String s){
        new AnnounceFrame(s);
        //throw new IllegalArgumentException(s);
    }
    @Override
    public void clean() {
        if (token != null) {
            token.stop();
            token.exit();
        }
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
        Nxy = FFTUtils.bestDimension((int)(sizeXY + paddingSizeXY.getValue()));
        Nz= FFTUtils.bestDimension((int)( sizeZ + paddingSizeZ.getValue()));
        outputShape = Shape.make(Nxy, Nxy, Nz);
    }

    private void setDefaultValue() {
        weightsMethod.setValue( weightOptions[3]);
        radial.setValue(false);
        if(debug){
            resultCostPrior.setValue(   "No results yet");
            resultDefocus.setValue(     "No results yet");
            resultModulus.setValue(     "No results yet");
            resultPhase.setValue(       "No results yet");
        }

        if (!isHeadLess()) {
            outputSize.setEnabled(false);
            imageSize.setEnabled(false);
            mu.setEnabled(false);
            if(debug){
                resultCostPrior.setEnabled(false);
                resultDefocus.setEnabled(false);
                resultModulus.setEnabled(false);
                resultPhase.setEnabled(false);
            }
        }
    }

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

        imageGlob = new EzPanel("Step 1: Data"); //Border layout to be sure that the images are stacked to the up
        EzPanel imagePan = new EzPanel("FILEPanel");
        image = new EzVarSequence("Sequence:");
        channel = new EzVarChannel("Canal:", image.getVariable(), false);
        imageSize = new EzVarText("Image size:");
        outputSize = new EzVarText("Output size:");
        paddingSizeXY = new EzVarInteger("padding xy:",30, Integer.MAX_VALUE,1);
        paddingSizeZ = new EzVarInteger("padding z :",30, Integer.MAX_VALUE,1);
        saveMetaData = new EzButton("Save metadata", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if (debug) {
                    System.out.println("Saving metadata");
                }
                updateMetaData();
            }
        });

        expertMode.addVarChangeListener(new EzVarListener<Boolean>() {
            @Override
            public void variableChanged(EzVar<Boolean> source, Boolean newValue) {
            	paddingSizeXY.setVisible(newValue);
            	paddingSizeZ.setVisible(newValue);
            	varianceGlob.setVisible(newValue);
            	epsilon.setVisible(newValue);
            	show.setVisible(newValue);
            	
            	nbAlphaCoef.setVisible(newValue); 
                nbBetaCoef.setVisible(newValue);   
                radial.setVisible(newValue);
                defocusMaxIter.setVisible(newValue);
                phaseMaxIter.setVisible(newValue);
                modulusMaxIter.setVisible(newValue); 
            }
        });


        updatePaddedSize();
        updateOutputSize();   
        updateImageSize();


        EzVarListener<Integer> zeroPadActionListener = new EzVarListener<Integer>() {
            @Override
            public void variableChanged(EzVar<Integer> source, Integer newValue) {
                updatePaddedSize();
                updateImageSize();
                updateOutputSize();  
                resetPSF();
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
                if (seq != null) {
                    meta = getMetaData(seq);
                    sizeX = seq.getSizeX();
                    sizeY = seq.getSizeY();
                    sizeZ = seq.getSizeZ();
                    dxy_nm.setValue(    meta.dxy);
                    dz_nm.setValue(     meta.dz);
                    na.setValue(     meta.na);
                    lambda.setValue( meta.lambda);
                    ni.setValue(     meta.ni);
                    updatePaddedSize();
                    updateOutputSize();
                    updateImageSize();

                    imageShape = Shape.make(sizeX, sizeY, sizeZ);
                    if (debug) {
                        System.out.println("Seq changed:" + sizeX + "  "+ Nxy);
                    }
                    // setting restart value to the current sequence
                    restart.setValue(newValue);
                }
            }

        });
        image.setNoSequenceSelection();


        EzVarListener<Double> metaActionListener = new EzVarListener<Double>() {
            @Override
            public void variableChanged(EzVar<Double> source, Double newValue) {
                Sequence seq = image.getValue();
                if (seq != null)  {              
                 //   setMetaData(seq) ;
                }
                resetPSF();
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
        lambda = new EzVarDouble( "\u03BB(nm):",542.,10.,15000.,10);
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

        // Note The listener of PSF is after BDEC tab

        /****************************************************/
        /**                WEIGHTING TAB                   **/
        /****************************************************/
        //Creation of the inside of WEIGHTING TAB
        varianceGlob = new EzPanel("Step 1b: Weights"); //Border layout to be sure that the images are stacked to the up
        EzPanel varianceTab = new EzPanel("VarianceTab");
        weightsMethod = new EzVarText(      "Weighting:", weightOptions, false);
        weights = new EzVarSequence(        "Map:");
        gain = new EzVarDouble(             "Gain:",1.,0.01,Double.MAX_VALUE,1);
        noise = new EzVarDouble(            "Readout Noise:",10.,0.,Double.MAX_VALUE,0.1);
        deadPixGiven = new EzVarBoolean(    "Bad data?", false);
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

        deadPixGiven.addVarChangeListener(new EzVarListener<Boolean>() {
            @Override
            public void variableChanged(EzVar<Boolean> source, Boolean newValue) {
                deadPixel.setVisible(deadPixGiven.getValue());
                deadPixel.setNoSequenceSelection();
            }
        });

        weights.setVisible(false);
        gain.setVisible(false);
        noise.setVisible(false);
        deadPixel.setVisible(false);
        showWeight = new EzButton("Show weight map", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                showWeightClicked();
                if (debug) {
                    System.out.println("Weight compute");
                }
            }
        });	

        
        /****************************************************/
        /**                    DECONV TAB                  **/
        /****************************************************/
        //Creation of the inside of DECONVOLUTION TAB
        deconvGlob = new EzPanel("Step 2: Deconvolution"); //Border layout to be sure that the images are stacked to the up
        EzPanel deconvTab = new EzPanel("DeconvolutionTab");
        mu = new EzVarDouble("Regularization level:",1E-5,0.,Double.MAX_VALUE,0.01);
        logmu = new EzVarDouble("Log10 of the Regularization level:",-5,-Double.MAX_VALUE,Double.MAX_VALUE,1);
        epsilon = new EzVarDouble("Threshold level:",1E-2,0.,Double.MAX_VALUE,0.01);
        nbIteration = new EzVarInteger("Number of iterations: ",10,0,Integer.MAX_VALUE ,1);
        positivity = new EzVarBoolean("Enforce nonnegativity:", true);
        // crop = new EzVarBoolean("Crop output to match input:", false);
        restart = new EzVarSequence("Starting point:");
        docDec = new EzLabel("Launch a deconvolution using the current PSF", Color.red);
        startDec = new EzButton("Start Deconvolution", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                Thread workerThread = new Thread() {
                    public void run() {
                        launch(true);
                    }
                };
                workerThread.start();
            }
        });
        restart.setNoSequenceSelection();
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
        
        cropResultInDeconv = new EzButton(          "Show deconvolved image at the size of the input", new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) { //FIXME Use View
                // We crop on the fly the last tvdec result to match input
                if (tvDec != null) {
                    ShapedArray prevResult = tvDec.getResult();
                    if (prevResult != null) {
                        Sequence croppedResult = new Sequence("Cropped Result");
                        ShapedArray croppedArray = ArrayUtils.crop(prevResult, imageShape);
                        double[] in = croppedArray.toDouble().flatten();
                        for (int j = 0; j < sizeZ; j++) {
                            double[] temp = new double[sizeX*sizeY];
                            for (int i = 0; i < sizeX*sizeY; i++) {
                                temp[i] = in[i+j*sizeX*sizeY];
                            }
                            croppedResult.setImage(0,j, new IcyBufferedImage(sizeX, sizeY, temp));
                        }
                        // TODO add meta data
                        setMetaData(croppedResult) ;
                        addSequence(croppedResult);
                    }
                }
            }
        });

        EzGroup groupStop1 = new EzGroup("Emergency STOP", stopDec);

        /****************************************************/
        /**                      BDEC TAB                  **/
        /****************************************************/
        //Creation of the inside of BDec TAB
        bdecGlob = new EzPanel("Step 3: Blind dec."); //Border layout to be sure that the images are stacked to the up
        EzPanel bdecTab = new EzPanel("BDecTab");
        nbAlphaCoef = new EzVarText(            "Number of phase coefs N\u03B1:", nAlphaOptions, 3,false );
        nbBetaCoef = new EzVarText(             "Number of modulus coefs N\u03B2:", nBetaOptions,0, false);       
        radial = new EzVarBoolean(              "Radially symmetric PSF", false);
        defocusMaxIter = new EzVarInteger(      "Max. nb. of iterations for defocus:",20,0,Integer.MAX_VALUE ,1);
        phaseMaxIter = new EzVarInteger(        "Max. nb. of iterations for phase:",20,0,Integer.MAX_VALUE ,1);
        modulusMaxIter = new EzVarInteger(      "Max. nb. of iterations for modulus:",0,0,Integer.MAX_VALUE ,1);
        bDecTotalIteration = new EzVarInteger(  "Number of loops:",2,0,Integer.MAX_VALUE ,1);

        resetPSF = new EzButton(            "Reset PSF", new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                resetPSF();
            }
        });

        radial.addVarChangeListener(new EzVarListener<Boolean>() {
            @Override
            public void variableChanged(EzVar<Boolean> source, Boolean newValue) {
                if(newValue){
                    nbAlphaCoef.setDefaultValues(         nAlphaOptionsR,3, false );
                    nbBetaCoef.setDefaultValues(         nBetaOptionsR,1, false );
                    resetPSF();
                }else{
                    nbAlphaCoef.setDefaultValues(         nAlphaOptions,7, false );
                    nbBetaCoef.setDefaultValues(         nBetaOptions,1, false );
                    resetPSF();           		
                }

            }
        });

        psfShow2 = new EzButton(        "Show PSF", new ActionListener() {
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
                    public void run() {
                        launch(false);
                    }
                };
                workerThread.start();
            }
        });

         show = new EzGroup("PSF visualization", psfShow2, showPhase, showModulus);
        stopBlind = new EzButton("STOP Computation", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                stopExecution();
            }
        });

        EzGroup groupStop2 = new EzGroup("Emergency STOP", stopBlind);

        cropResult = new EzButton(          "Show deconvolved image at the size of the input", new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                // We crop on the fly the last tvdec result to match input
                if (tvDec != null) {
                    ShapedArray prevResult = tvDec.getResult();
                    if (prevResult != null) {
                        Sequence croppedResult = new Sequence("Cropped Result");
                        ShapedArray croppedArray = ArrayUtils.crop(prevResult, imageShape);
                        double[] in = croppedArray.toDouble().flatten();
                        for (int j = 0; j < sizeZ; j++) {
                            double[] temp = new double[sizeX*sizeY];
                            for (int i = 0; i < sizeX*sizeY; i++) {
                                temp[i] = in[i+j*sizeX*sizeY];
                            }
                            croppedResult.setImage(0,j, new IcyBufferedImage(sizeX, sizeY, temp));
                        }
                        // TODO add meta data
                        setMetaData(croppedResult) ;
                        addSequence(croppedResult);
                    }
                }
            }
        });


        /****************************************************/
        /**                    RESULT TAB                  **/
        /****************************************************/
        resultGlob = new EzPanel("Results"); //Border layout to be sure that the images are stacked to the up
        EzPanel resultTab = new EzPanel(    "ResultsTab");
        resultCostPrior = new EzVarText(    "Costs");
        resultDefocus = new EzVarText(      "Defocus");
        resultModulus = new EzVarText(      "Modulus");
        resultPhase = new EzVarText(        "Phase");



        /****************************************************/
        /**                      ToolTips                  **/
        /****************************************************/
        image.setToolTipText(ToolTipText.sequenceImage);
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
        nbIteration.setToolTipText(ToolTipText.doubleMaxIter);
        paddingSizeXY.setToolTipText(ToolTipText.doublePadding);
        paddingSizeZ.setToolTipText(ToolTipText.doublePadding);

        nbAlphaCoef.setToolTipText(ToolTipText.doubleNalpha);
        nbBetaCoef.setToolTipText(ToolTipText.doubleNbeta);
        //grtolPhase.setToolTipText(ToolTipText.doubleGrtolPhase);
        //grtolModulus.setToolTipText(ToolTipText.doubleGrtolModulus);
        //grtolDefocus.setToolTipText(ToolTipText.doubleGrtolDefocus);
        bDecTotalIteration.setToolTipText(ToolTipText.doubleBDecTotalIteration);

        restart.setToolTipText(ToolTipText.booleanRestart);
        positivity.setToolTipText(ToolTipText.booleanPositivity);
        cropResultInDeconv.setToolTipText(ToolTipText.booleanCrop);
        cropResult.setToolTipText(ToolTipText.booleanCrop);
        resultTab.setToolTipText(ToolTipText.textOutput);

        /******** Adding ********/
        paddingSizeXY.setVisible(false);
    	paddingSizeZ.setVisible(false);
    	varianceGlob.setVisible(false);
    	epsilon.setVisible(false);
    	
    	nbAlphaCoef.setVisible(false); 
        nbBetaCoef.setVisible(false);   
        radial.setVisible(false);
        defocusMaxIter.setVisible(false);
        phaseMaxIter.setVisible(false);
        modulusMaxIter.setVisible(false); 
        
        /**** IMAGE ****/
        imagePan.add(image);
        imagePan.add(channel);
        imagePan.add(imageSize);
        imagePan.add(paddingSizeXY);
        imagePan.add(paddingSizeZ);
        imagePan.add(outputSize);
        imagePan.add(dxy_nm);
        imagePan.add(dz_nm);


        imagePan.add(na);
        imagePan.add(ni);
        imagePan.add(lambda);     

        imagePan.add(saveMetaData);
        imagePan.add(showPSF);

        imageGlob.add(imagePan);
        tabbedPane.add(imageGlob);

        /**** Variance ****/
        varianceTab.add(weightsMethod);
        varianceTab.add(weights);
        varianceTab.add(gain);
        varianceTab.add(noise);
        varianceTab.add(deadPixGiven);
        varianceTab.add(deadPixel);
        varianceTab.add(showWeight);
        //Creation of VARIANCE TAB
        varianceGlob.add(varianceTab);
        tabbedPane.add(varianceGlob);

        /**** Deconv ****/

        deconvTab.add(logmu);
        deconvTab.add(mu);
        deconvTab.add(epsilon);
        deconvTab.add(nbIteration);
        deconvTab.add(positivity);
        deconvTab.add(restart);
        deconvTab.add(docDec);
        deconvTab.add(startDec);
        deconvTab.add(cropResultInDeconv);
        deconvTab.add(groupStop1);
        //Creation of DECONVOLUTION TAB
        deconvGlob.add(deconvTab);
        tabbedPane.add(deconvGlob);

        /**** BDec ****/
        bdecTab.add(resetPSF);
        bdecTab.add(nbAlphaCoef);
        bdecTab.add(nbBetaCoef);
        bdecTab.add(radial);
        bdecTab.add(defocusMaxIter);
        bdecTab.add(phaseMaxIter);
        bdecTab.add(modulusMaxIter);
        bdecTab.add(bDecTotalIteration);
        bdecTab.add(docBlind);
        bdecTab.add(startBlind);
        bdecTab.add(cropResult);
        bdecTab.add(show);
        bdecTab.add(groupStop2);
        //Creation of BDec TAB
        bdecGlob.add(bdecTab);
        tabbedPane.add(bdecGlob);

        if(debug){
            /**** Result ****/
            resultTab.add(resultCostPrior);
            resultTab.add(resultDefocus);
            resultTab.add(resultModulus);
            resultTab.add(resultPhase);
            resultGlob.add(resultTab);
            tabbedPane.add(resultGlob);
        }

        addEzComponent(expertMode);
        addEzComponent(tabbedPane);
        // Must be added to global panel first 
        show.setFoldedState(true);

        setDefaultValue();
        if (isHeadLess()) {
            outputHeadlessImage = new EzVarSequence("Output Image");
            outputHeadlessPSF = new EzVarSequence("Output PSF");
        }

        token = new ReconstructionThreadToken(new double[]{mu.getValue(),epsilon.getValue(),0.0,grtol});
        thread = new ReconstructionThread(token);
        thread.start();
    }

    protected void resetPSF() {
        defocuSpace = null;
        defocusVector= null;
        alphaSpace = null;
        alphaVector = null;
        betaSpace = null;
        betaVector = null;
        buildpupil();
    }
    public boolean launchDeconvolution(DoubleArray imgArray, DoubleArray psfArray, DoubleArray weight){
        return launchDeconvolution(imgArray, psfArray, weight, true, false);
    }

    public boolean launchDeconvolution(DoubleArray imgArray, DoubleArray psfArray, DoubleArray weight, boolean cleanPrevResult, boolean ignoreRestart){
        if (tvDec == null) {
            tvDec = new TotalVariationJobForIcy(token);
            tvDec.setResult(null);
        }
        tvDec.setAbsoluteTolerance(0.0);
        tvDec.setWeight(weight);
        tvDec.setData(imgArray);
        tvDec.setPsf(psfArray);
        tvDec.setViewer(new TvViewer());
        thread.setJob(tvDec);
        // If We are at step 0 and a previous result was given, we give the result
        // Else if explicitly asked we clean 
        // Else we use previous result
        if (!ignoreRestart && restart.getValue() != null) {
            Sequence restartSeq = restart.getValue();
            // We verify that the previous result is conform to our expectations: !Null and same dim as input
            if (restartSeq != null) {
                int numCanal = channel.getValue();
                ShapedArray tmpDoubleArray = IcyBufferedImageUtils.imageToArray(restartSeq, numCanal);
                boolean sameAsOrigin = tmpDoubleArray.getRank() == 3 && tmpDoubleArray.getDimension(0) == sizeX 
                        && tmpDoubleArray.getDimension(1) == sizeY && tmpDoubleArray.getDimension(2) == sizeZ;
                boolean sameAsPrevious = tmpDoubleArray.getRank() == 3 && tmpDoubleArray.getDimension(0) == Nxy 
                        && tmpDoubleArray.getDimension(1) == Nxy && tmpDoubleArray.getDimension(2) == Nz;
                if (!(sameAsOrigin || sameAsPrevious)) {
                    throwError("The previous result does not have the same dimensions as the input image");
                    return false;
                }
                tmpDoubleArray = ArrayUtils.pad(tmpDoubleArray, outputShape);
                tvDec.setResult(tmpDoubleArray);
            }
        } else if(cleanPrevResult){
            tvDec.setResult(null);
        } else {
            tvDec.setResult(tvDec.getResult());
        }
        tvDec.setPositivity(positivity.getValue());
        tvDec.setRegularizationWeight(mu.getValue());
        tvDec.setRegularizationThreshold(epsilon.getValue());
        tvDec.setRelativeTolerance(grtol);
        tvDec.setMaximumIterations(nbIteration.getValue());
        tvDec.setOutputShape(outputShape);
        token.start();
        setResult(tvDec);
        return true;
    }

    @Override
    protected void execute() {
        launch(false);
    }

    protected void launch(boolean runDeconv) {
        try {
            startBlind.setText("Computing...");
            if (isHeadLess()) { // For a trigger to update all values
                image.valueChanged(null, null, null);
            }
            if (debug) {
                System.out.println("-------------IMAGE-------------------");
                System.out.println("File: "+image.getValue());              //Used
                System.out.println("Canal: "+channel.getValue());
                System.out.println("--------------PSF------------------");
                //         System.out.println("PSF: "+psf.getValue());                 //Used
                System.out.println("dxy: "+dxy_nm.getValue()*1E-9);
                System.out.println("dz: "+dz_nm.getValue()*1E-9);
                System.out.println("Nxy: "+Nxy);
                System.out.println("Nx: "+Nz);
                System.out.println("NA: "+na.getValue());
                System.out.println("\u03BB: "+lambda.getValue()*1E-9);
                System.out.println("ni: "+ni.getValue());
                System.out.println("--------------Variance------------------");
                System.out.println("Weights: "+weights.getValue());
                System.out.println("Gain: "+gain.getValue());
                System.out.println("Noise: "+noise.getValue());
                System.out.println("deadPix: "+deadPixel.getValue());
                System.out.println("--------------DECONV------------------");
                System.out.println("zeroPad xy: "+paddingSizeXY.getValue());
                System.out.println("zeroPad z: "+paddingSizeZ.getValue());
                System.out.println("nbIter: "+nbIteration.getValue());
                System.out.println("Restart: "+restart.getValue());
                System.out.println("Positivity: "+positivity.getValue());
                System.out.println("--------------BDEC------------------");
                System.out.println("nbIter: "+nbIteration.getValue());
                System.out.println("zeroPad: "+paddingSizeXY.getValue());
                /*System.out.println("nbIterZern: "+grtolPhase.getValue());
                System.out.println("module: "+grtolModulus.getValue());
                System.out.println("defoc: "+grtolDefocus.getValue());*/
                System.out.println("Number of total iterations: "+bDecTotalIteration.getValue());
                System.out.println("------------------------------------");
                System.out.println("");
            }
            run = true;

            // Preparing parameters and testing input
            Sequence imgSeq = image.getValue();
            if (imgSeq == null)
            {
                throwError("An image/sequence of images should be given");
                return;
            }

            // Set the informations about the input
            if (sizeZ == 1) {
                throwError("Input data must be 3D");
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
            runBdec = !runDeconv;


            if (pupil == null) {
                buildpupil();
            }

            imgArray = (DoubleArray) IcyBufferedImageUtils.imageToArray(imgSeq, imageShape, numCanal);
            wgtArray = createWeights(imgArray).toDouble();
            
            // Zero pad the data and their weights. FIXME: should not be needed (to save memory).
            imgArray = (DoubleArray) ArrayUtils.pad(imgArray, outputShape);
            wgtArray = (DoubleArray) ArrayUtils.pad(wgtArray, outputShape);



            /*---------------------------------------*/
            /*            OPTIMISATION               */
            /*---------------------------------------*/

            if (runBdec) {
                if ((alphaVector==null)||( Integer.parseInt(nbAlphaCoef.getValue()) != nbAlpha)){
                    nbAlpha = Integer.parseInt(nbAlphaCoef.getValue());
                    if (nbAlpha==0){
                        guessPhase = false; 	
                    }else{
                        guessPhase = true;
                        alphaSpace = new DoubleShapedVectorSpace(new int[]{nbAlpha});
                        alphaVector = alphaSpace.create();
                    }
                }
                if  ((betaVector==null)||(Integer.parseInt(nbBetaCoef.getValue()) != nbBeta)){
                    nbBeta = Integer.parseInt(nbBetaCoef.getValue());
                    if (nbBeta==0){
                        guessModulus = false; 	
                    }else{
                        guessModulus = true;
                        double[] beta = new double[nbBeta];
                        beta[0] = 1;
                        betaSpace = new DoubleShapedVectorSpace(new int[]{beta.length});
                        betaVector = betaSpace.wrap(beta);
                    }
                }

                if (defocuSpace==null){
                    double[] defocus = {ni.getValue()/(lambda.getValue()*1E-9), 0., 0.}; 
                    defocuSpace = new DoubleShapedVectorSpace(new int[]{defocus.length});
                    defocusVector = defocuSpace.wrap(defocus);
                }

                PSFEstimation = new PSF_Estimation(pupil);

                PSFEstimation.setWeight(wgtArray);
                PSFEstimation.setData(imgArray);

                PSFEstimation.enablePositivity(false);
                PSFEstimation.setAbsoluteTolerance(0.0);

                for(int i = 0; i < bDecTotalIteration.getValue(); i++) {
                     psfArray = (DoubleArray) ArrayUtils.roll(Double3D.wrap(pupil.getPSF(), outputShape));
                     pupil.freePSF();
                    /* OBJET ESTIMATION (by the current PSF) */
                    // If first iteration we use given result, after we continue with our previous result (i == 0)
                    if (!launchDeconvolution(imgArray, psfArray, wgtArray, false, !(i == 0))) {
                        return;
                    }
                    PSFEstimation.setObj(tvDec.getResult());

                    /* Defocus estimation */
                    if (defocusMaxIter.getValue()>0){
                        if (debug && verbose) {
                            System.out.println("------------------");
                            System.out.println("Defocus estimation");
                            System.out.println("------------------");
                        }
                        PSFEstimation.setRelativeTolerance(0.);
                        PSFEstimation.setMaximumIterations(defocusMaxIter.getValue());
                        PSFEstimation.fitPSF(defocusVector, PSF_Estimation.DEFOCUS);
                        System.out.println( "Defocus   "+Arrays.toString(PSFEstimation.getPupil().getDefocusMultiplyByLambda())  );
                    }

                    /* Phase estimation */
                    if((phaseMaxIter.getValue()>0)& guessPhase){
                        if (debug && verbose) {
                            System.out.println("Phase estimation");
                            System.out.println("------------------");
                        }
                        PSFEstimation.setResult(null);                    
                        PSFEstimation.setMaximumIterations(phaseMaxIter.getValue());
                        PSFEstimation.fitPSF(alphaVector, PSF_Estimation.ALPHA);
                    }


                    /* Modulus estimation */
                    if((modulusMaxIter.getValue() >0)&guessModulus){
                        if (debug && verbose) {
                            System.out.println("Modulus estimation");
                            System.out.println("------------------");
                        }
                        PSFEstimation.setResult(null);           
                        PSFEstimation.setMaximumIterations(modulusMaxIter.getValue());
                        PSFEstimation.fitPSF(betaVector, PSF_Estimation.BETA);
                        // MathUtils.normalise(betaVector.getData());
                    }
                    if (debug || isHeadLess()) {
                        showResult(i);
                    }

                    //If we want a emergency stop
                    if (!run) {
                        return;
                    }
                }
                pupil = PSFEstimation.getPupil();
            } else {
                 psfArray = (DoubleArray) ArrayUtils.roll(Double3D.wrap(pupil.getPSF(), outputShape));
                 pupil.freePSF();
                launchDeconvolution(imgArray, psfArray, wgtArray);
            }
            pupil.freePSF();
            // Everything went well, the restart will be the current sequence
            // For now, if not called in the graphic thread launch errors
            // I'am using lastsequence because invokelater will find null with sequence
            SwingUtilities.invokeLater(new Runnable() {
                public void run() {
                    restart.setValue(lastSequence);
                }
            });
            // In any cases the next image will be in a new sequence
            lastSequence = sequence;
            sequence = null;
        } catch (IllegalArgumentException e) {
            new AnnounceFrame("Oops, Error: "+ e.getMessage());
            if (debug) {
                e.printStackTrace();
            }
        } finally {
            startBlind.setText("Guess PSF");
        }
    }

    /**
     * The goal is to create an array of weights, but it will be created depending
     * the user input so we will have to test each cases:
     *  	- None
     *  	- A given map
     *  	- A variance map
     *  	- A computed variance
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
        
        if (deadPixGiven.getValue() && (seq = deadPixel.getValue()) != null) {
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

    /**
     * Get the results from a reconstruction and plot the intermediate result
     * 
     * @param tvDec
     */
    private void setResult(ReconstructionJob tvDec){
        try{
            //Here we will update the sequence
            if (sequence == null || (sequence != null && sequence.isEmpty())) {
                sequence = new Sequence();
                setMetaData(image.getValue(), sequence);
                if (isHeadLess()) {
                    outputHeadlessImage.setValue(sequence);
                } else {
                    addSequence(sequence);
                }
            }
            sequence.beginUpdate();
            double[] in = tvDec.getResult().toDouble().flatten();
            for (int j = 0; j < Nz; j++) {
                double[] temp = new double[Nxy*Nxy];
                for (int i = 0; i < Nxy*Nxy; i++) {
                    temp[i] = in[i+j*Nxy*Nxy];
                }
                sequence.setImage(0,j, new IcyBufferedImage(Nxy, Nxy, temp));
            }


            sequence.endUpdate();
            sequence.setName("TV mu:"+mu.getValue()+" Iteration:"+tvDec.getIterations());

            System.out.println("Cost "+tvDec.getCost() );
            //Then we will update the result tab panel
            if (runBdec) {
                if (debug) {
                    resultCostPrior.setValue("Cost "+tvDec.getCost()                       );
                    resultDefocus.setValue(   "Defocus   "+Arrays.toString(pupil.getDefocusMultiplyByLambda())   );
                    resultModulus.setValue(  "Modulus    "+pupil.getRho()[0]                     );
                    resultPhase.setValue(    "Phase      "+pupil.getPhi()[0]                     );
                }
            }
        } catch (NullPointerException e) {
            //Here in case of brutal stop the sequence can become null but it's not important as it's an emergency stop
            //So we do nothing
            System.out.println("INFO: Emergency stop detected in setResult");
            e.printStackTrace();
        }
    }

    /**
     * A debug function linked to the PSF
     * 
     * @param num Number of PSF already generated
     */
    private void showResult(int num)
    {
        Sequence psf3DSequence = new Sequence();
        DoubleArray psf = Double3D.wrap(pupil.getPSF(), outputShape);
        double[] PSF_shift = ArrayUtils.roll(psf).toDouble().flatten();
        //double[] PSF_shift = MathUtils.fftShift3D(pupil.getPSF(), xyPad, xyPad, sizeZPad);
        psf3DSequence.setName("PSF Estimated - " + num);
        for (int k = 0; k < Nz; k++)
        {
            psf3DSequence.setImage(0, k, new IcyBufferedImage(Nxy, Nxy,
                    MathUtils.getArray(PSF_shift, Nxy, Nxy, k)));
        }
        if (isHeadLess()) {
            outputHeadlessPSF.setValue(psf3DSequence);
        } else {
            addSequence(psf3DSequence);
        }

    }

    /*****************************************/
    /** All the PSF buttons call are here   **/
    /*****************************************/

    private void buildpupil()
    {
        pupil = new WideFieldModel(na.getValue(), lambda.getValue()*1E-9, ni.getValue(), dxy_nm.getValue()*1E-9,
                dz_nm.getValue()*1E-9, Nxy, Nxy, Nz,radial.getValue());
    }	

    private void psfClicked()
    {
    	/* PSF0 initialisation */
    	if(pupil==null)
    	{
    		buildpupil();
    	}

    	/* PSF0 Sequence */
    	Sequence PSF0Sequence = new Sequence();

    	DoubleArray psf = Double3D.wrap(pupil.getPSF(), outputShape);
    	double[] PSF_shift = ArrayUtils.roll(psf).toDouble().flatten();
    	//double[] PSF_shift = MathUtils.fftShift3D(pupil.getPSF(), xyPad, xyPad, sizeZPad);
    	for (int k = 0; k < Nz; k++)
    	{
    		PSF0Sequence.setImage(0, k, new IcyBufferedImage(Nxy, Nxy,
    				MathUtils.getArray(PSF_shift, Nxy, Nxy, k)));
    	}
    	setMetaData(PSF0Sequence) ;

    	PSF0Sequence.setName("PSF");
    	addSequence(PSF0Sequence);
    }

    
    private void phaseClicked()
    {
        /* PSF0 initialisation */ 
        if(pupil==null)
        {
            buildpupil();
        }
        /* Phase Sequence */
        Sequence phaseSequence = new Sequence();
        phaseSequence.setName("Phase of the pupil");
        DoubleArray psf = Double2D.wrap(pupil.getPhi(), Shape.make(Nxy, Nxy));
        double[] phase_shift = ArrayUtils.roll(psf).toDouble().flatten();
        //double[] phase_shift = MathUtils.fftShift1D(pupil.getPhi(), xyPad, xyPad);
        phaseSequence.addImage(new IcyBufferedImage(Nxy, Nxy, phase_shift));
        addSequence(phaseSequence);
    }

    private void modulusClicked()
    {
        /* PSF0 initialisation */ 
        if(pupil==null)
        {
            buildpupil();
        }
        /* Modulus Sequence */
        Sequence modulusSequence = new Sequence();
        modulusSequence.setName("Modulus of the pupil");
        DoubleArray psf = Double2D.wrap(pupil.getRho(), Shape.make(Nxy, Nxy));
        double[] modulus_shift = ArrayUtils.roll(psf).toDouble().flatten();
        //double[] modulus_shift = MathUtils.fftShift1D(pupil.getRho(), xyPad, xyPad);
        modulusSequence.addImage(new IcyBufferedImage(Nxy, Nxy, modulus_shift));
        addSequence(modulusSequence);
    }

    private void showWeightClicked()
    {

        //    DoubleArray weight = createWeights(...).toDouble();
        /* PSF0 Sequence */
        //  FIXME fix the weightsize
        Sequence img = image.getValue();
        if (img == null) {
            new AnnounceFrame("No image input found");
            return;
        }
        Shape myShape=Shape.make(sizeX, sizeY, sizeZ);

        int numCanal = channel.getValue();
        DoubleArray input = (DoubleArray) IcyBufferedImageUtils.imageToArray(img, myShape, numCanal);
        Sequence WeightSequence = new Sequence();
        WeightSequence.setName("Weight");
        double[] wght = createWeights(input).toDouble().flatten();
        if (wght.length != sizeX*sizeY*sizeZ) {
            new AnnounceFrame("Invalid weight size");
            return;
        }
        for (int k = 0; k < sizeZ; k++)
        {
            WeightSequence.setImage(0, k, new IcyBufferedImage(sizeX, sizeY,
                    MathUtils.getArray(wght, sizeX, sizeY, k)));
        }
        addSequence(WeightSequence);
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
    private MyMetadata getMetaData(Sequence seq){
        OMEXMLMetadata metDat = seq.getMetadata();
        if (meta == null) {
            meta = new MyMetadata();
            if (metDat.getInstrumentCount() > 0) {
                try {
                    meta.na      = metDat.getObjectiveLensNA(0, 0);
                    //meta.lambda  = metDat.getChannelEmissionWavelength(0, 0).getValue().doubleValue()*1E6;  //I suppose the value I will get is in um
                } catch(Exception e){
                    System.out.println("Failed to get some metadatas, will use default values for na, lambda");
                }
            } else {
                if (debug && verbose) {
                    System.out.println("INFO: Metadata: No instrument so no metadata.");
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

    //Copy the input metadata to the output. We may want to change some with our values
    //So it should be done here
    private void setMetaData(Sequence seqOld, Sequence seqNew) {
        OMEXMLMetadataImpl newMetdat = OMEUtil.createOMEMetadata(seqOld.getMetadata());
        //newMetdat.setImageDescription("MyDescription", 0);
        newMetdat.setPixelsPhysicalSizeX(OMEUtil.getLength(dxy_nm.getValue()*1E-9), 0);
        newMetdat.setPixelsPhysicalSizeY(OMEUtil.getLength(dxy_nm.getValue()*1E-9), 0);
        newMetdat.setPixelsPhysicalSizeZ(OMEUtil.getLength(dz_nm.getValue()*1E-9), 0);
        // seqNew.setMetaData(newMetdat); //FIXME may not working now
    }


    //Copy the input metadata to the output. We may want to change some with our values
    //So it should be done here
    private void setMetaData(Sequence seqNew) {
        OMEXMLMetadataImpl newMetdat = OMEUtil.createOMEMetadata();
        //newMetdat.setImageDescription("MyDescription", 0);
        newMetdat.setPixelsPhysicalSizeX(OMEUtil.getLength(dxy_nm.getValue()*1E-3), 0);
        newMetdat.setPixelsPhysicalSizeY(OMEUtil.getLength(dxy_nm.getValue()*1E-3), 0);
        newMetdat.setPixelsPhysicalSizeZ(OMEUtil.getLength(dz_nm.getValue()*1E-3), 0);
        seqNew.setMetaData((OMEXMLMetadataImpl) newMetdat); //FIXME may not working now
    }

    private void updateMetaData() {
        Sequence seq = image.getValue();
        if (seq != null) {
            try {
                OMEXMLMetadata newMetdat = MetaDataUtil.generateMetaData(seq, false);
                newMetdat.setPixelsPhysicalSizeX(OMEUtil.getLength(dxy_nm.getValue()*1E-3), 0);
                newMetdat.setPixelsPhysicalSizeY(OMEUtil.getLength(dxy_nm.getValue()*1E-3), 0);
                newMetdat.setPixelsPhysicalSizeZ(OMEUtil.getLength(dz_nm.getValue()*1E-3), 0);
                seq.setMetaData((OMEXMLMetadataImpl) newMetdat); //FIXME may not working now
            } catch (ServiceException e) {
                e.printStackTrace();
            }
        } else {
            new AnnounceFrame("Nothing to save");
        }
    }

    @Override
    public void stopExecution() {
        if (token != null) {
            token.stop();
        }
        if (PSFEstimation != null) {
            PSFEstimation.stop();
        }
        run = false;
    }

    //Debug function Will have to be deleted in the future 
    @SuppressWarnings("unused")
    private void ZernikeClicked()
    {
        /* PSF0 initialisation */
        if(pupil==null)
        {
            buildpupil();
        }

        /* PSF0 Sequence */
        Sequence PSF0Sequence = new Sequence();
        	 
            for (int k = 0; k < pupil.getNZern(); k++)
            {
                PSF0Sequence.setImage(0, k, new IcyBufferedImage(Nxy, Nxy,
                       ArrayUtils.roll(Double2D.wrap(pupil.getZernike(k), Shape.make(Nxy,Nxy))).toDouble().flatten()));
            }
        PSF0Sequence.setName("Zernike");
        addSequence(PSF0Sequence);
    }
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
    @Override
    public void declareInput(VarList inputMap) {// FIXME Use subclasses for protocols
        initialize();
        inputMap.add("image", image.getVariable());
        inputMap.add("channel", channel.getVariable());
        inputMap.add("dxy", dxy_nm.getVariable());
        inputMap.add("dz", dz_nm.getVariable());
        inputMap.add("NA", na.getVariable());
        inputMap.add("ni", ni.getVariable());
        inputMap.add("lambda", lambda.getVariable());
        
        inputMap.add("deadPixel", deadPixel.getVariable());
        inputMap.add("gain", gain.getVariable());
        inputMap.add("noise", noise.getVariable());
        
        inputMap.add("mu", mu.getVariable());
        inputMap.add("espilon", epsilon.getVariable());
        inputMap.add("Postivity", positivity.getVariable());
        inputMap.add("nbIteration", nbIteration.getVariable());
        inputMap.add("restart", restart.getVariable());
        
        inputMap.add("radial",radial.getVariable());
        inputMap.add("nbAlphaCoef", nbAlphaCoef.getVariable());
        inputMap.add("nbBetaCoef", nbBetaCoef.getVariable());
        inputMap.add("defocusMaxIter", defocusMaxIter.getVariable());
        inputMap.add("phaseMaxIter", phaseMaxIter.getVariable());
        inputMap.add("modulusMaxIter", modulusMaxIter.getVariable());
        inputMap.add("bDecTotalIteration", bDecTotalIteration.getVariable());

    }
    @Override
    public void declareOutput(VarList outputMap) {
    	outputMap.add("outputSize", outputSize.getVariable());
        outputMap.add("outputImage", outputHeadlessImage.getVariable());
        outputMap.add("outputPSF", outputHeadlessPSF.getVariable());
    }
	@Override
	public String getMainPluginClassName() {
		
		// TODO Auto-generated method stub
		return MitivDeconvolution.class.getName();
	}
}