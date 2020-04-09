/*
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

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.SwingUtilities;

import icy.file.Loader;
import icy.gui.frame.progress.AnnounceFrame;
import icy.main.Icy;
import icy.plugin.PluginDescriptor;
import icy.plugin.PluginInstaller;
import icy.plugin.PluginRepositoryLoader;
import icy.plugin.PluginUpdater;
import icy.sequence.Sequence;
import icy.system.thread.ThreadUtil;
import icy.util.StringUtil;
import mitiv.array.ArrayUtils;
import mitiv.array.DoubleArray;
import mitiv.array.ShapedArray;
import mitiv.base.Shape;
import mitiv.base.Traits;
import mitiv.conv.WeightedConvolutionCost;
import mitiv.cost.DifferentiableCostFunction;
import mitiv.cost.HyperbolicTotalVariation;
import mitiv.jobs.DeconvolutionJob;
import mitiv.linalg.shaped.DoubleShapedVectorSpace;
import mitiv.linalg.shaped.FloatShapedVectorSpace;
import mitiv.linalg.shaped.ShapedVectorSpace;
import mitiv.utils.FFTUtils;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzButton;
import plugins.adufour.ezplug.EzGroup;
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
import plugins.mitiv.io.DeconvHook;
import plugins.mitiv.io.IcyImager;


/**
 * SimpleDEMIC implements a TV regularized multi-dimensional deconvolution.
 * @author Ferréol
 *
 */
public class SimpleDEMIC extends DEMICSPlug implements Block, EzStoppable {


    /***************************************************/
    /**         Plugin interface variables            **/
    /***************************************************/

    private EzVarSequence psf; // Point Spread Function

    private EzVarInteger    paddingSizeX,paddingSizeY;

    /** deconvolution tab: **/
    private EzVarDouble  epsilon;
    private EzVarBoolean  showIteration;
    private EzVarBoolean  normalizePSF;

    private EzButton saveParam, loadParam;
    /** headless mode: **/



    /****************************************************/
    /**                 VARIABLES                      **/
    /****************************************************/

    static boolean debug =false;

    private int psfSizeX=1,psfSizeY=1,psfSizeZ=1;

    private EzVarChannel channelpsf;
    private EzGroup ezPaddingGroup;
    private EzGroup ezWeightingGroup;
    private EzGroup ezDeconvolutionGroup;
    private EzGroup ezDeconvolutionGroup2;

    private DeconvolutionJob deconvolver;

    private ShapedArray badArray=null;

    private ShapedVectorSpace dataSpace, objectSpace;

    /*********************************/
    /**      Initialization         **/
    /*********************************/
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

        data = new EzVarSequence("Data:");
        channel = new EzVarChannel("Data channel:", data.getVariable(), false);
        psf = new EzVarSequence("PSF:");
        channelpsf = new EzVarChannel("PSF channel :", psf.getVariable(), false);
        normalizePSF = new EzVarBoolean("Normalize PSF (sum=1):", true);
        restart = new EzVarSequence("Starting point:");
        restart.setNoSequenceSelection();
        channelRestart = new EzVarChannel("Initialization channel :", restart.getVariable(), false);
        dataSize = new EzVarText("Data size:");
        outputSize = new EzVarText("Output size:");
        paddingSizeX = new EzVarInteger("padding x:",0, Integer.MAX_VALUE,1);
        paddingSizeY = new EzVarInteger("padding y:",0, Integer.MAX_VALUE,1);
        paddingSizeZ = new EzVarInteger("padding z :",0, Integer.MAX_VALUE,1);
        dx_nm = new EzVarDouble("dx(nm):",64.5,0., Double.MAX_VALUE,1.);
        dy_nm = new EzVarDouble("dy(nm):",64.5,0., Double.MAX_VALUE,1.);
        dz_nm = new EzVarDouble("dz(nm):",64.5,0., Double.MAX_VALUE,1.);

        dataSize.setVisible(false);
        outputSize.setVisible(false);

        data.setNoSequenceSelection();
        psf.setNoSequenceSelection();
        ezPaddingGroup = new EzGroup("Padding",paddingSizeX,paddingSizeY,paddingSizeZ);
        ezPaddingGroup.setFoldedState(true);

        EzVarListener<Integer> zeroPadActionListener = new EzVarListener<Integer>() {
            @Override
            public void variableChanged(EzVar<Integer> source, Integer newValue) {
                updateImageSize();
                updatePaddedSize();
            }
        };
        paddingSizeX.addVarChangeListener(zeroPadActionListener);
        paddingSizeY.addVarChangeListener(zeroPadActionListener);
        paddingSizeZ.addVarChangeListener(zeroPadActionListener);




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
        data.addVarChangeListener(new EzVarListener<Sequence>() {
            @Override
            public void variableChanged(EzVar<Sequence> source,
                    Sequence newValue) {
                if (debug) {
                    System.out.println("Seq ch..."+data.getValue());
                }
                dataChanged() ;
                updatePSFSize();

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
        weightsMethod = new EzVarText(      "Weighting:", weightOptions,3, false);
        weights = new EzVarSequence(        "Map:");
        gain = new EzVarDouble(             "Gain:",1.,0.01,Double.MAX_VALUE,1);
        noise = new EzVarDouble(            "Readout Noise:",10.,0.,Double.MAX_VALUE,0.1);
        deadPixel = new EzVarSequence(      "Bad data map:");
        weights.setNoSequenceSelection();

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
                    weights.setVisible(false);
                    gain.setVisible(false);
                    noise.setVisible(false);
                } else if (StringUtil.equals(weightsMethod.getValue() , weightOptions[1] )|| StringUtil.equals(weightsMethod.getValue() , weightOptions[2])) {  //Personnalized map ou Variance map
                    weights.setVisible(true);
                    gain.setVisible(false);
                    noise.setVisible(false);
                    weights.setNoSequenceSelection();
                } else if (StringUtil.equals(weightsMethod.getValue() , weightOptions[3])) {  //Computed variance
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
                DoubleArray dataArray;
                // Preparing parameters and testing input

                Sequence dataSeq = data.getValue();
                if(dataSeq!=null){
                    dataArray =    sequenceToArray(dataSeq, channel.getValue()).toDouble();
                    wgtArray = createWeights(dataArray,badArray).toDouble();
                    IcyImager.show(wgtArray,null,"Weight map",false);
                }
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
        epsilon = new EzVarDouble("Threshold level:",1E-2,0.,Double.MAX_VALUE,1.0);
        nbIterDeconv = new EzVarInteger("Number of iterations: ",10,1,Integer.MAX_VALUE ,1);
        positivity = new EzVarBoolean("Enforce nonnegativity:", true);
        singlePrecision = new EzVarBoolean("Compute in single precision:", false);
        showIteration = new EzVarBoolean("Show intermediate results:", true);

        if (isHeadLess()){
            showIteration.setValue(false);
        }

        scale = new EzVarDoubleArrayNative("Aspect ratio of a voxel",  new double[][] { {1.0 ,1.0, 1.0} },true);
        ezDeconvolutionGroup2 = new EzGroup("More  parameters",epsilon,scale,positivity,singlePrecision);
        ezDeconvolutionGroup2.setFoldedState(true);
        ezDeconvolutionGroup = new EzGroup("Deconvolution parameters",logmu,mu,nbIterDeconv,ezDeconvolutionGroup2);

        startDec = new EzButton("Start Deconvolution", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if ( deconvolver!=null && deconvolver.isRunning()){
                    stopExecution();
                }else{
                    Thread workerThread = new Thread() {
                        @Override
                        public void run() {
                            launch();


                        }
                    };
                    workerThread.start();
                }
            }
        });

        showFullObject = new EzButton("Show the full (padded) object", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(debug){
                    System.out.println("showFullObject");
                }
                Sequence fSeq;
                fSeq = new Sequence("Deconvolved image");
                fSeq.copyMetaDataFrom(data.getValue(), false);
                //     show(solver.getSolution(),fSeq,"Deconvolved "+ data.getValue().getName() + " mu="+solver.getRegularizationLevel() );
                if(objArray != null){
                    IcyImager.show(objArray,fSeq,"Deconvolved "+ data.getValue().getName() + "with padding. mu="+ mu.getValue() ,isHeadLess());
                }else {
                    IcyImager.show(ArrayUtils.extract(dataArray, outputShape),fSeq,"Deconvolved "+ data.getValue().getName() + "with padding. mu="+ mu.getValue(),isHeadLess() );
                }
            }
        });

        EzVarListener<Double> logmuActionListener = new EzVarListener<Double>() {
            @Override
            public void variableChanged(EzVar<Double> source, Double newValue) {
                mu.setValue(Math.pow(10, logmu.getValue()));
            }
        };
        logmu.addVarChangeListener(logmuActionListener);

        EzGroup groupVisu  = new EzGroup("Visualization", showIteration,showFullObject);
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


        loadParam = new EzButton("Load parameters", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(saveFile.getValue()!=null){
                    loadParameters(saveFile.getValue());
                    Sequence dataSeq = data.getValue();
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
                    }                }
                if (debug) {
                    System.out.println("Load parameters");
                }
            }
        });


        addEzComponent(data);
        addEzComponent(channel);
        addEzComponent(psf);
        addEzComponent(channelpsf);
        addEzComponent(normalizePSF);
        addEzComponent(restart);
        addEzComponent(channelRestart);
        addEzComponent(dataSize);
        addEzComponent(outputSize);
        addEzComponent(ezPaddingGroup);
        addEzComponent(ezWeightingGroup);
        addEzComponent(ezDeconvolutionGroup);
        addEzComponent(startDec);
        addEzComponent(groupVisu);
        addEzComponent(saveFile);
        addEzComponent(saveParam);
        addEzComponent(loadParam);

        setDefaultValue();

        updatePaddedSize();
        updateOutputSize();
        updateImageSize();

        if (!isHeadLess()) {
            outputSize.setEnabled(false);
            dataSize.setEnabled(false);
            mu.setEnabled(false);
        }else{
            outputHeadlessImage = new EzVarSequence("Output Image");
            outputHeadlessWght = new EzVarSequence("Computed weight");
        }


    }



    /**
     *
     */
    private void saveParamClicked() {
        File pathName = saveFile.getValue();
        if(pathName!=null){
            if (!pathName.getName().endsWith(".xml")){
                pathName = new File(pathName.getAbsolutePath()+".xml");
            }
            saveParameters(pathName);
        }
    }




    @Override
    protected void updatePaddedSize() {

        Nx = FFTUtils.bestDimension(Math.max(psfSizeX, sizeX + paddingSizeX.getValue()));
        Ny = FFTUtils.bestDimension(Math.max(psfSizeY,sizeY + paddingSizeY.getValue()));
        if ((Nz==1)&&(paddingSizeZ.getValue()==0)){
            outputShape = new Shape(Nx, Ny);
        }else{
            Nz= FFTUtils.bestDimension(Math.max(psfSizeZ,sizeZ + paddingSizeZ.getValue()));
            outputShape = new Shape(Nx, Ny, Nz);
        }
        updateOutputSize();
        if(debug){
            System.out.println(" UpdatePaddedSize" + paddingSizeX.getValue()  + outputShape.toString());
        }
    }

    protected void updatePSFSize() {
        paddingSizeX.setValue(Math.max(psfSizeX/10,10));
        paddingSizeY.setValue(Math.max(psfSizeY/10,10));
        if( (psfSizeZ==1))
            paddingSizeZ.setValue(  0);
        else
            paddingSizeZ.setValue(  Math.max(psfSizeZ/10,10));

        updatePaddedSize();
        if(debug){
            System.out.println(" UpdatePaddedSize " + paddingSizeX.getValue()  + outputShape.toString());
        }
    }

    /****************************************************/
    /**                  RUN PLUGIN                    **/
    /****************************************************/
    @Override
    protected void execute() {

        if (isHeadLess()){
            if(  Icy.getCommandLinePluginArgs().length!=0){
                initialize();
                parseCmdLine();
            }
            showIteration.setValue(false);
            if(debug){
                System.out.println("Launch it:"+nbIterDeconv.getValue());
            }
        }
        launch();

    }

    protected void launch() {
        //       solver = new EdgePreservingDeconvolution();

        try {
            startDec.setText("Emergency stop");

            dataSeq = data.getValue();
            Sequence psfSeq = psf.getValue();
            Sequence restartSeq = restart.getValue();


            if (dataSeq == null)
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
            if ((sizeZ == 1)&&(sizeY == 1)) {
                throwError("Input data must be 2D or 3D");
                return;
            }
            if (paddingSizeX.getValue() < 0.0) {
                throwError("Padding value cannot be negative");
                return;
            }
            if (paddingSizeY.getValue() < 0.0) {
                throwError("Padding value cannot be negative");
                return;
            }
            if (paddingSizeZ.getValue() < 0.0) {
                throwError("Padding value cannot be negative");
                return;
            }

            dataArray =  sequenceToArray(dataSeq, channel.getValue());
            psfArray =  sequenceToArray(psfSeq,  channelpsf.getValue());
            dataShape = dataArray.getShape();
            if (restart.getValue() != null && restartSeq != null){
                objArray =  sequenceToArray(restartSeq, channelRestart.getValue());
                if(debug){
                    System.out.println("restart seq:" +restartSeq.getName());
                }
            }else{
                objArray = sequenceToArray(dataSeq, channel.getValue());
                if(debug){
                    System.out.println("restart seq is null:");
                }
            }

            if(singlePrecision.getValue()){
                wgtArray = createWeights(dataArray.toFloat(),badArray).toFloat();
            }else{
                wgtArray = createWeights(dataArray.toDouble(),badArray).toDouble();
            }
            cursequence = new Sequence("Current Iterate");
            cursequence.copyMetaDataFrom(dataSeq, false);


            if(Nz==1){
                if (scale.getValue().length !=2){
                    throwError("Pixel scale must have 2 elements");
                    return;
                }
            }else{
                if (scale.getValue().length !=3){
                    throwError("Pixel scale must have 3 elements");
                    return;
                }
            }
            if (debug){
                System.out.println("Launch it:"+nbIterDeconv.getValue());
            }

            IcyImager curImager = new IcyImager(cursequence, isHeadLess());

            DeconvHook dHook = new DeconvHook(curImager, dataShape,null, debug);
            DeconvHook dHookfinal = new DeconvHook(curImager, dataShape,"Deconvolved "+dataSeq.getName(), debug);
            // deconvolver = new EdgePreservingDeconvolutionJob(dataArray, psfArray, wgtArray, outputShape, mu.getValue(), epsilon.getValue(), scale.getValue(), positivity.getValue(), singlePrecision.getValue(), nbIterDeconv.getValue(), dHook , dHookfinal);

            buildVectorSpaces();

            DifferentiableCostFunction fprior = new HyperbolicTotalVariation(objectSpace, epsilon.getValue(), scale.getValue());
            WeightedConvolutionCost fdata =  WeightedConvolutionCost.build( objectSpace, dataSpace);
            fdata.setData(dataArray);
            fdata.setWeights(wgtArray);
            fdata.setPSF(psfArray, normalizePSF.getValue());
            deconvolver  = new DeconvolutionJob( fdata,  mu.getValue(),fprior,  positivity.getValue(),nbIterDeconv.getValue(),  dHook,  dHookfinal);

            objArray = ArrayUtils.extract(objArray, outputShape, 0.); //Padding to the right size
            objArray = deconvolver.deconv(objArray);


            SwingUtilities.invokeLater(new Runnable() {
                @Override
                public void run() {
                    if(debug){
                        System.out.println("invoke later");
                    }
                    restart.setValue(cursequence);
                    channelRestart.setValue(0);

                    if (isHeadLess()) {
                        if(outputHeadlessImage==null){
                            outputHeadlessImage  = new EzVarSequence("Output Image");
                        }
                        outputHeadlessImage.setValue(cursequence);
                        if (outputHeadlessWght==null) {
                            outputHeadlessWght = new EzVarSequence("Computed weights");
                        }
                        outputHeadlessWght.setValue(arrayToSequence(wgtArray));

                        if(outputPath!=null){
                            IcyImager.save(cursequence, outputPath);
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
            startDec.setText("Start Deconvolution");
        }
    }





    /**
     *
     */
    private void buildVectorSpaces() {
        /* Determine the floating-point type for all vectors. */
        int type;
        if (singlePrecision.getValue()) {
            type = Traits.FLOAT;
        } else if (dataArray.getType() == Traits.DOUBLE ||
                (psfArray != null && psfArray.getType() == Traits.DOUBLE) ||
                (wgtArray != null && wgtArray.getType() == Traits.DOUBLE)) {
            type = Traits.DOUBLE;
        } else {
            type = Traits.FLOAT;
        }
        /* Build vector spaces. */
        if (type == Traits.FLOAT) {
            dataSpace = new FloatShapedVectorSpace(dataShape);
            objectSpace = new FloatShapedVectorSpace(outputShape);
        } else {
            dataSpace = new DoubleShapedVectorSpace(dataShape);
            objectSpace = new DoubleShapedVectorSpace(outputShape);
        }
    }



    //If the user call the stop button
    @Override
    public void stopExecution() {
        if (deconvolver !=null)
            deconvolver.abort();
    }
    /**
     *  set default values of the plugin
     */
    @Override
    protected void setDefaultValue() {
        super.setDefaultValue();
        paddingSizeX.setValue(10);
        paddingSizeY.setValue(10);
        paddingSizeZ.setValue(10);
    }

    //The input variable for the protocol
    @Override
    public void declareInput(VarList inputMap) {
        initialize();

        super.declareInput(inputMap);
        inputMap.add("psf", psf.getVariable());
        inputMap.add("psf channel", channelpsf.getVariable());

        inputMap.add("Padding X", paddingSizeX.getVariable());
        inputMap.add("Padding Y", paddingSizeY.getVariable());
        inputMap.add("Padding Z", paddingSizeZ.getVariable());


        inputMap.add("espilon", epsilon.getVariable());
    }



    @Override
    public void clean() {
        stopExecution();

    }

    private void parseCmdLine(){
        String[] args = Icy.getCommandLinePluginArgs();

        loadParameters( new File(args[0]));
        for (int i = 1; i < args.length; i++) {
            switch (args[i]) {
                case "-i":
                    if(i+1 >= args.length)
                        break;

                    System.out.println("load image:" + args[i+1]);
                    data.setValue(Loader.loadSequence(args[i+1], 0, false));
                    data.valueChanged(data.getVariable(), null, data.getValue());
                    if(i+3 >= args.length)
                        break;
                    if( args[i+1].startsWith("-"))
                        break;
                    if(args[i+2].equalsIgnoreCase("-c")){
                        channel.setValue(Integer.parseInt(args[i+3]));
                        i=i+3;
                    }else{
                        i++;
                    }


                    break;
                case "-p":
                    if(i+1 >= args.length)
                        break;
                    if( args[i+1].startsWith("-"))
                        break;
                    System.out.println("load psf:" + args[i+1]);
                    psf.setValue(Loader.loadSequence(args[i+1], 0, false));
                    if(i+3 >= args.length)
                        break;
                    if(args[i+2].equalsIgnoreCase("-c")){
                        channelpsf.setValue(Integer.parseInt(args[i+3]));
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
                    outputPath = args[i+1];
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
                    System.out.println("-p psf file");
                    System.out.println("-r restart file");
                    System.out.println("-o deconvolved output file");
                    System.out.println("-badpix bad pixels file");
                    System.out.println("-wghtmap weight or variance map file");
                    break;
            }
        }


    }


}