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

import static plugins.mitiv.io.Icy2TiPi.sequenceToArray;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import icy.gui.frame.progress.AnnounceFrame;
import icy.gui.frame.progress.FailedAnnounceFrame;
import icy.plugin.PluginDescriptor;
import icy.plugin.PluginInstaller;
import icy.plugin.PluginRepositoryLoader;
import icy.plugin.PluginUpdater;
import icy.sequence.MetaDataUtil;
import icy.sequence.Sequence;
import icy.system.thread.ThreadUtil;
import icy.util.OMEUtil;
import icy.util.StringUtil;
import mitiv.array.ArrayFactory;
import mitiv.array.ArrayUtils;
import mitiv.array.ByteArray;
import mitiv.array.ShapedArray;
import mitiv.base.Shape;
import mitiv.base.Traits;
import mitiv.conv.WeightedConvolutionCost;
import mitiv.cost.DifferentiableCostFunction;
import mitiv.jobs.DeconvolutionJob;
import mitiv.linalg.shaped.DoubleShapedVectorSpace;
import mitiv.linalg.shaped.FloatShapedVectorSpace;
import mitiv.linalg.shaped.ShapedVectorSpace;
import mitiv.utils.FFTUtils;
import mitiv.utils.HistoMap;
import mitiv.weights.WeightFactory;
import mitiv.weights.WeightUpdater;
import ome.units.UNITS;
import ome.xml.meta.OMEXMLMetadata;
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
import plugins.adufour.ezplug.EzVarFile;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarListener;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.ezplug.EzVarText;
import plugins.mitiv.io.IcyImager;

/**
 * Plugin class for all plugins of DEconvolution MIcroscopy Studio
 * (DEMICS)
 * @author Ferréol
 *
 */
public abstract class DEMICSPlug extends EzPlug  implements Block,EzStoppable{

    protected EzVarSequence   dataEV;           // data
    protected EzVarChannel    channelEV;        // data channel

    protected  OMEXMLMetadata metDat = null;

    protected EzVarDouble     logmu, mu;      // deconvolution hyper parameters; mu = 10^(logmu)
    protected EzVarSequence   restartEV;        // starting point
    protected EzVarChannel    channelRestartEV; // starting point channel
    protected EzVarBoolean    positivityEV;     // enforce non negativity
    protected EzButton        startDecButton, showFullObjectButton;


    protected EzVarText       dataSizeTxt;       //
    protected EzVarText       outputSizeTxt;     // size of the object after padding
    // optical parameters
    protected EzVarDouble     dy_nm,dx_nm, dz_nm;  //  pixels size in (x,y) and z
    protected EzVarDouble     na=null  ;             //  numerical aperture
    protected EzVarDouble     lambda=null;         //  wavelength
    protected EzVarDouble     ni=null;             //  refractive index of the immersion index

    protected EzVarInteger    nbIterDeconv;   // number of iteration for the deconvolution stage
    protected EzVarBoolean    singlePrecision;// compute in single precision
    protected EzVarDoubleArrayNative scale;   // scale of a voxel should be [1 1 dz/dxy]

    protected EzVarInteger    paddingSizeXY, paddingSizeZ; // number of pixels added in each direction

    // Main variables for the deconvolution part
    protected EzGroup ezDeconvolutionGroup;
    protected int sizeX=128, sizeY=128, sizeZ=64; // Input sequence sizes
    protected  int Nx=128,Ny=128, Nz=64;             // Output (padded sequence size)
    protected Shape psfShape = new Shape(Nx, Ny, Nz);
    protected Shape outputShape;
    protected Sequence dataSeq=null;
    protected Sequence cursequence; // Sequence containing the current solution
    protected Shape dataShape;
    protected ShapedArray wgtArray, dataArray, psfArray, objArray;
    protected ShapedArray modelArray=null;
    protected ByteArray badpixArray=null;

    protected EzGroup ezWeightingGroup, groupVisu;
    protected EzVarText       weightsMethod;  // Combobox for variance estimation
    protected final String[] weightOptions = new String[]{"Variance =1","Inverse variance map","Variance map","Computed variance","Automatic variance estimation"};
    protected EzVarDouble     gain, noise;    // gain of the detector in e-/lvl and detector noise in e-
    protected EzVarSequence weightsSeq, badpixMap; // maps of inverse variance and bad pixels
    protected EzButton        showWeightButton;

    protected EzGroup ezPaddingGroup;


    protected EzVarFile       saveFile, loadFile;// xml files to save and load parameters
    protected EzVarBoolean    showIteration;  // show object update at each iteration
    protected EzVarSequence   outputHeadlessImage=null;
    protected EzVarSequence   outputHeadlessWght=null;
    protected EzButton saveParam, loadParam;
    protected String outputPath=null;

    protected DeconvolutionJob deconvolver ;
    protected ShapedVectorSpace dataSpace, objectSpace;
    protected  int vectorSpaceType;
    protected DifferentiableCostFunction fprior;
    protected WeightedConvolutionCost fdata;

    // listener
    protected EzVarListener<Integer> zeroPadActionListener;
    protected WeightUpdater wghtUpdt=null;
    /*********************************/
    /**            DEBUG            **/
    /*********************************/
    private boolean debug = false;      // Show psf steps
    @SuppressWarnings("unused")
    private boolean verbose = false;    // Show some values, need debug to true


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


        dataEV = new EzVarSequence("Data:");
        channelEV = new EzVarChannel("Data channel:", dataEV.getVariable(), false);
        restartEV = new EzVarSequence("Starting point:");
        restartEV.setNoSequenceSelection();
        channelRestartEV = new EzVarChannel("Initialization channel :", restartEV.getVariable(), false);
        dataSizeTxt = new EzVarText("Data size:");
        outputSizeTxt = new EzVarText("Output size:");
        paddingSizeZ = new EzVarInteger("padding in z :",0, Integer.MAX_VALUE,1);
        dx_nm = new EzVarDouble("dx(nm):",64.5,0., Double.MAX_VALUE,1.);
        dy_nm = new EzVarDouble("dy(nm):",64.5,0., Double.MAX_VALUE,1.);
        dz_nm = new EzVarDouble("dz(nm):",64.5,0., Double.MAX_VALUE,1.);

        dataSizeTxt.setVisible(false);
        outputSizeTxt.setVisible(false);
        dataEV.setNoSequenceSelection();

        zeroPadActionListener = new EzVarListener<Integer>() {
            @Override
            public void variableChanged(EzVar<Integer> source, Integer newValue) {
                updateImageSize();
                updatePaddedSize();
            }
        };
        paddingSizeZ.addVarChangeListener(zeroPadActionListener);




        restartEV.addVarChangeListener(new EzVarListener<Sequence>() {
            @Override
            public void variableChanged(EzVar<Sequence> source,
                    Sequence newValue) {
                newValue = restartEV.getValue();
                if(debug){
                    if (newValue != null || (newValue != null && newValue.isEmpty())) {
                        System.out.println("restart changed:"+newValue.getName());
                    }
                }
            }
        });

        channelEV.addVarChangeListener(new EzVarListener<Integer>() {
            @Override
            public void variableChanged(EzVar<Integer> source, Integer newValue) {
                channelRestartEV.setValue(newValue);
                dataChanged();
                modelArray = null;
            }
        });




        /****************************************************/
        /**                WEIGHTING GROUP                   **/
        /****************************************************/
        weightsMethod = new EzVarText(      "Weighting:", weightOptions,4, false);
        weightsSeq = new EzVarSequence(        "Map:");
        gain = new EzVarDouble(             "Gain:",1.,Double.MIN_VALUE,Double.MAX_VALUE,0.1);
        noise = new EzVarDouble(            "Readout Noise:",10.,0.,Double.MAX_VALUE,0.1);
        badpixMap = new EzVarSequence(      "Bad data map:");
        weightsSeq.setNoSequenceSelection();

        weightsMethod.addVarChangeListener(new EzVarListener<String>() {
            @Override
            public void variableChanged(EzVar<String> source, String newValue) {
                if(debug){
                    System.out.println("weight:" + weightsMethod.getValue()+".");
                    System.out.println("weight:" + newValue+".");
                    System.out.println("weight:" + weightOptions[3]+".");
                    System.out.println("weight:" + weightOptions[3]==newValue);
                }
                if (StringUtil.equals(weightsMethod.getValue(), weightOptions[0])) { //None
                    weightsSeq.setVisible(false);
                    gain.setVisible(false);
                    noise.setVisible(false);
                } else if (StringUtil.equals(weightsMethod.getValue() , weightOptions[1] )
                        || StringUtil.equals(weightsMethod.getValue() , weightOptions[2])) {  //Personnalized map ou Variance map
                    weightsSeq.setVisible(true);
                    gain.setVisible(false);
                    noise.setVisible(false);
                    weightsSeq.setNoSequenceSelection();
                } else if (StringUtil.equals(weightsMethod.getValue() , weightOptions[3])
                        || StringUtil.equals(weightsMethod.getValue() , weightOptions[4])) {  //Computed variance
                    weightsSeq.setVisible(false);
                    gain.setVisible(true);
                    noise.setVisible(true);
                    weightsSeq.setNoSequenceSelection();
                } else {
                    throwError("Invalid argument passed to weight method");
                    return;
                }
            }
        });


        weightsSeq.setVisible(false);
        gain.setVisible(false);
        noise.setVisible(false);
        badpixMap.setVisible(true);
        badpixMap.setNoSequenceSelection();
        showWeightButton = new EzButton("Show weight map", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                // Preparing parameters and testing input
                Sequence dataSeq = dataEV.getValue();
                if(wgtArray!=null) {
                    IcyImager.show(wgtArray,null,"Weight map",false);
                }else if(dataSeq!=null){
                    dataArray =    sequenceToArray(dataSeq, channelEV.getValue()).toDouble();
                    createWeights(true);
                    IcyImager.show(wgtArray,null,"Weight map",false);
                }
                if (debug) {
                    System.out.println("Weight compute");

                }
            }
        });




        /****************************************************/
        /**                    DECONV GROUP                  **/
        /****************************************************/
        mu = new EzVarDouble("Regularization level:",1,0.,Double.MAX_VALUE,0.01);
        logmu = new EzVarDouble("Log10 of the Regularization level:",0,-Double.MAX_VALUE,Double.MAX_VALUE,1);
        nbIterDeconv = new EzVarInteger("Number of iterations: ",50,1,Integer.MAX_VALUE ,1);
        positivityEV = new EzVarBoolean("Enforce nonnegativity:", true);
        singlePrecision = new EzVarBoolean("Compute in single precision:", false);
        showIteration = new EzVarBoolean("Show intermediate results:", true);
        scale = new EzVarDoubleArrayNative("Aspect ratio of a voxel",  new double[][] { {1.0 ,1.0, 1.0} },true);

        if (isHeadLess()){
            showIteration.setValue(false);
        }



        showFullObjectButton = new EzButton("Show the full (padded) object", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(debug){
                    System.out.println("showFullObjectButton");
                }
                Sequence fSeq;
                fSeq = new Sequence("Deconvolved image");
                fSeq.copyMetaDataFrom(dataEV.getValue(), false);
                if(objArray != null){
                    IcyImager.show(objArray,fSeq,"Deconvolved "+ dataEV.getValue().getName() + "with padding. mu="+ mu.getValue() ,isHeadLess());
                }else {
                    IcyImager.show(ArrayUtils.extract(dataArray, outputShape),fSeq,"data "+ dataEV.getValue().getName() + "with padding. mu="+ mu.getValue(),isHeadLess() );
                }
            }
        });

        EzVarListener<Double> logmuActionListener = new EzVarListener<Double>() {
            @Override
            public void variableChanged(EzVar<Double> source, Double newValue) {
                mu.setValue(Math.pow(10, logmu.getValue()));
                mu.setEnabled(false);
            }
        };
        logmu.addVarChangeListener(logmuActionListener);

        groupVisu  = new EzGroup("Visualization", showIteration,showFullObjectButton);
        groupVisu.setFoldedState(true);

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

    }
    /**
     *
     */
    protected abstract void loadParamClicked();
    /**
     *
     */
    protected void saveParamClicked() {
        File pathName = saveFile.getValue();
        if(pathName!=null){
            if (!pathName.getName().endsWith(".xml")){
                pathName = new File(pathName.getAbsolutePath()+".xml");
            }
            saveParameters(pathName);
        }
    }

    /**
     * The goal is to create an array of weights, but it will be created depending
     * the user input so we will have to test each cases:
     *      - None
     *      - A given map
     *      - A variance map
     *      - A computed variance
     * Then we apply the dead pixel map
     */

    protected void createWeights( boolean normalize) {

        badpixArray = (ByteArray) ArrayFactory.create(Traits.BYTE, dataShape);
        if (WeightFactory.flagBads(dataArray,badpixArray,dataSeq.getChannelTypeMax(channelEV.getValue()))){
            if (!isHeadLess()) {
                new AnnounceFrame("Warning, saturated pixel detected, accounting them as dead pixels", "show", new Runnable() {
                    @Override
                    public void run() {
                        Sequence deadSequence = new Sequence("Saturations map");
                        deadSequence.copyMetaDataFrom(dataSeq, false);
                        IcyImager.show(badpixArray, deadSequence, "saturations map", isHeadLess());
                    }
                }, 3);
            }
        }
        if (weightsMethod.getValue() == weightOptions[1]) {
            // A map of weights is provided.
            Sequence seq;
            if ((seq = weightsSeq.getValue()) != null) {
                wgtArray =  sequenceToArray(seq);
                if (!wgtArray.equals(dataArray)){
                    throwError("Weight map must have the same size than the data");
                }
            }
        } else if (weightsMethod.getValue() == weightOptions[2]) {
            // A variance map is provided. FIXME: check shape and values.
            Sequence seq;
            if ((seq = weightsSeq.getValue()) != null) {
                ShapedArray varArray =  sequenceToArray(seq);
                if (!varArray.equals(dataArray)){
                    throwError("Variance map must have the same size than the data");
                }
                wgtArray = WeightFactory.computeWeightsFromVariance(varArray);
            }
        } else if (weightsMethod.getValue() == weightOptions[3]) {
            // Weights are computed given the gain and the readout noise of the detector.
            double gamma = gain.getValue();
            double sigma = noise.getValue();
            double alpha = gamma;
            double beta = (sigma/gamma)*(sigma/gamma);
            wgtArray = WeightFactory.computeWeightsFromData(dataArray,   alpha,  beta);
        } else if (weightsMethod.getValue() == weightOptions[4]) {
            // Weights are computed from the current estimate of the data modelArray =object * PSF
            //  the gain and the readout noise of the detector are automatically estimated from the variance of the data given the modelArray
            if (modelArray==null) { // without modelArray (before any deconvolution) rely on the"compute variance" method
                double gamma = gain.getValue();
                double sigma = noise.getValue();
                double alpha = gamma;
                double beta = (sigma/gamma)*(sigma/gamma);
                wgtArray = WeightFactory.computeWeightsFromData(dataArray,   alpha,  beta);
            }else {
                // wgtArray = WeightFactory.computeWeightsFromModel(dataArray,modelArray,badpixArray);
                HistoMap hm = new HistoMap(modelArray, dataArray, badpixArray);
                gain.setValue(hm.getAlpha());
                noise.setValue(Math.sqrt(hm.getBeta())/hm.getAlpha());
                wgtArray = hm.computeWeightMap(modelArray);
            }
        }

        WeightFactory.removeBads(wgtArray,  badpixArray);

        if (normalize) {
            WeightFactory.normalize(wgtArray);
        }
    }



    protected void createWeights() {
        createWeights(  false);
    }


    /**
     * Function triggered when the data sequence change.
     * Change parameters according to metadata of the data sequence
     *
     */
    protected void dataChanged() {
        dataSizeTxt.setVisible(false);
        outputSizeTxt.setVisible(false);
        dataSeq = dataEV.getValue();
        if(dataSeq!=null){
            sizeX = dataSeq.getSizeX();
            sizeY = dataSeq.getSizeY();
            sizeZ = dataSeq.getSizeZ();
            // setting restart value to the current sequence
            restartEV.setValue(dataSeq);
            channelRestartEV.setValue(channelEV.getValue());
            metDat = dataSeq.getOMEXMLMetadata();
            dx_nm.setValue( dataSeq.getPixelSizeX()*1E3);
            dy_nm.setValue( dataSeq.getPixelSizeY()*1E3);
            dz_nm.setValue( dataSeq.getPixelSizeZ()*1E3);

            if (lambda != null) {
                try {
                    lambda.setValue(metDat.getChannelEmissionWavelength(0, channelEV.getValue()).value(UNITS.NANOMETER)
                            .doubleValue());
                } catch (Exception e) {
                    System.out.println("Failed to get wavelength from metadata, will use default values ");
                    lambda.setValue(500.0);
                }
            }
            if (na != null) {
                try {
                    na.setValue(metDat.getObjectiveLensNA(0, 0));
                } catch (Exception e) {
                    System.out.println("Failed to get numerical aperture from metadata, will use default values ");
                    na.setValue(1.4);
                }
            }
            if (ni != null) {
                try {
                    if (metDat.getObjectiveSettingsRefractiveIndex(0) != null)
                        ni.setValue(metDat.getObjectiveSettingsRefractiveIndex(0));
                    else {
                        System.out.println("Failed to get refractive index from metadata, will use default values ");
                        ni.setValue(1.518);
                    }
                } catch (Exception e) {
                    System.out.println("Failed to get refractive index from metadata, will use default values ");
                    ni.setValue(1.518);
                }
            }

            if (sizeZ==1) {
                dataShape = new Shape(sizeX, sizeY);
                scale.setValue(new double[]{1.0 , dx_nm.getValue()/ dy_nm.getValue()} );
            }else {
                dataShape = new Shape(sizeX, sizeY, sizeZ);
                scale.setValue(new double[]{1.0 , dx_nm.getValue()/ dy_nm.getValue(), dx_nm.getValue()/ dz_nm.getValue() } );
            }

            updatePaddedSize();
            updateOutputSize();
            updateImageSize();
            dataSizeTxt.setVisible(true);
            outputSizeTxt.setVisible(true);
            outputSizeTxt.setEnabled(false);
            dataSizeTxt.setEnabled(false);

            startDecButton.setEnabled(true);
        }else {
            startDecButton.setEnabled(false);
        }
    }



    /* (non-Javadoc)
     * The input variable for the protocol
     * @see plugins.adufour.blocks.lang.Block#declareInput(plugins.adufour.blocks.util.VarList)
     */
    @Override
    public void declareInput(VarList inputMap) {
        inputMap.add("image", dataEV.getVariable());
        inputMap.add("image channel", channelEV.getVariable());
        inputMap.add("starting point", restartEV.getVariable());
        channelRestartEV = new EzVarChannel("Initialization channel :", restartEV.getVariable(), false);

        inputMap.add("starting point channel", channelRestartEV.getVariable());

        inputMap.add("weights Method",weightsMethod.getVariable());
        inputMap.add("badpixMap", badpixMap.getVariable());
        inputMap.add("gain", gain.getVariable());
        inputMap.add("noise", noise.getVariable());

        inputMap.add("mu", mu.getVariable());
        inputMap.add("scale", scale.getVariable());

        inputMap.add("nbIteration", nbIterDeconv.getVariable());
        inputMap.add("positivity", positivityEV.getVariable());
        inputMap.add("single precision", singlePrecision.getVariable());

        saveFile = new EzVarFile("Save parameters in", "");
        inputMap.add("saveFile",  saveFile.getVariable());

    }

    /* (non-Javadoc)
     * output variable for the protocol
     * @see plugins.adufour.blocks.lang.Block#declareOutput(plugins.adufour.blocks.util.VarList)
     */
    @Override
    public void declareOutput(VarList outputMap) {
        outputMap.add("outputSizeTxt", outputSizeTxt.getVariable());
        outputMap.add("output", outputHeadlessImage.getVariable());
        outputMap.add("weightmap", outputHeadlessWght.getVariable());
    }



    /**
     *  set default values of the plugin
     */
    protected void setDefaultValue() {
        weightsMethod.setValue( weightOptions[4]);
        dataEV.setNoSequenceSelection();
        badpixMap.setNoSequenceSelection();


        paddingSizeZ.setValue(30);
        badpixMap.setNoSequenceSelection();

    }

    /**
     * Update the text indicating the size of the data
     */
    protected void updateImageSize() {
        String text ;
        if (Nz==1)
            text= sizeX+"x"+sizeY;
        else
            text= sizeX+"x"+sizeY+"x"+sizeZ;
        dataSizeTxt.setValue(text);
    }

    /**
     * update metaData of the data sequence using the value indicated in the plugin
     */
    protected void updateMetaData() {
        Sequence seq = dataEV.getValue();
        if (seq != null) {
            OMEXMLMetadata newMetdat = MetaDataUtil.generateMetaData(seq, true);
            newMetdat.setPixelsPhysicalSizeX(OMEUtil.getLength(dx_nm.getValue()*1E-3), 0);
            newMetdat.setPixelsPhysicalSizeY(OMEUtil.getLength(dy_nm.getValue()*1E-3), 0);
            newMetdat.setPixelsPhysicalSizeZ(OMEUtil.getLength(dz_nm.getValue()*1E-3), 0);
            newMetdat.setObjectiveLensNA(na.getValue(),0, 0);
            //    for (int i = 0; i <  newMetdat.getImageCount(); i++) {
            newMetdat.setChannelEmissionWavelength(OMEUtil.getLength(lambda.getValue()*1E-3),0, channelEV.getValue());
            newMetdat.setObjectiveSettingsRefractiveIndex(ni.getValue(), 0);
            //     }
            seq.setMetaData(newMetdat); //FIXME may not working now
        } else {
            new AnnounceFrame("Nothing to save");
        }
    }


    /**
     * Print the size of the deconvolved image  in the plugin
     * throwError if the number of pixel is to high to indexed with int
     */
    protected void updateOutputSize() {

        String text;
        if (Nz==1){
            text= Nx+"x"+Ny;
        }else{
            text= Nx+"x"+Ny+"x"+Nz;
        }
        outputSizeTxt.setValue(text);
        if((1.0*Nx*Ny*Nz)>Math.pow(2, 30)){
            throwError("Padded image is too large (>2^30)");
        }
    }


    /**
     * Update the size of the output image according to padding size
     * indicated in the plugging. This output size is rounded to the
     * next integer suitable for fast FFT computation
     *
     */
    protected void updatePaddedSize() {

        if (paddingSizeXY.getValue() < 0.0) {
            throwError("Padding value cannot be negative");
            return;
        }
        if (paddingSizeZ.getValue() < 0.0) {
            throwError("Padding value cannot be negative");
            return;
        }

        Nx = FFTUtils.bestDimension(sizeX + paddingSizeXY.getValue());
        Ny = FFTUtils.bestDimension(sizeY + paddingSizeXY.getValue());
        Nz= FFTUtils.bestDimension(sizeZ + paddingSizeZ.getValue());
        outputShape = new Shape(Nx, Ny, Nz);
        if(debug){
            System.out.println(" UpdatePaddedSize" + paddingSizeXY.getValue()  + outputShape.toString());
        }
    }

    /*
     * print error message
     *
     * @param s
     * the error message
     */
    protected void throwError(final String s) {
        if(isHeadLess()){
            throw new IllegalArgumentException(s);
        } else {
            new FailedAnnounceFrame(s);

        }

    }

    /**
     *
     */
    protected void buildVectorSpaces() {
        /* Determine the floating-point type for all vectors. */
        if (singlePrecision.getValue()) {
            vectorSpaceType = Traits.FLOAT;
        } else if (dataArray.getType() == Traits.DOUBLE ||
                (psfArray != null && psfArray.getType() == Traits.DOUBLE) ||
                (wgtArray != null && wgtArray.getType() == Traits.DOUBLE)) {
            vectorSpaceType = Traits.DOUBLE;
        } else {
            vectorSpaceType = Traits.FLOAT;
        }
        /* Build vector spaces. */
        if (vectorSpaceType == Traits.FLOAT) {
            dataSpace = new FloatShapedVectorSpace(dataShape);
            objectSpace = new FloatShapedVectorSpace(outputShape);
        } else {
            dataSpace = new DoubleShapedVectorSpace(dataShape);
            objectSpace = new DoubleShapedVectorSpace(outputShape);
        }
    }


}
