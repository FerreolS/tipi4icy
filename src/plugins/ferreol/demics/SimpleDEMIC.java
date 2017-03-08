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

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.SwingUtilities;

import icy.file.Loader;
import icy.main.Icy;
import icy.sequence.Sequence;
import mitiv.array.ArrayUtils;
import mitiv.array.DoubleArray;
import mitiv.base.Shape;
import mitiv.cost.EdgePreservingDeconvolution;
import mitiv.optim.OptimTask;
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


/**
 * SimpleDEMIC implements a TV regularized multi-dimensional deconvolution.
 * @author ferreol
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

    private EzButton saveParam, loadParam;
    /** headless mode: **/



    /****************************************************/
    /**                 VARIABLES                      **/
    /****************************************************/

    static boolean debug =false;

    private int psfSizeX=1,psfSizeY=1,psfSizeZ=1;
    private static double[][] scaleDef =new double[][] {{1.0},{1.0 ,1.0},{1.0 ,1.0, 1.0},{1.0 ,1.0, 1.0,1.0}};

    private EdgePreservingDeconvolution solver =  new EdgePreservingDeconvolution();
    private boolean run;
    private EzVarChannel channelpsf;
    private EzGroup ezPaddingGroup;
    private EzGroup ezWeightingGroup;
    private EzGroup ezDeconvolutionGroup;
    private EzGroup ezDeconvolutionGroup2;



    /*********************************/
    /**      Initialization         **/
    /*********************************/
    @Override
    protected void initialize() {

        if (!isHeadLess()) {
            getUI().setParametersIOVisible(false);
            getUI().setActionPanelVisible(false);
        }

        data = new EzVarSequence("Data:");
        channel = new EzVarChannel("Data channel:", data.getVariable(), false);
        psf = new EzVarSequence("PSF:");
        channelpsf = new EzVarChannel("PSF channel :", psf.getVariable(), false);
        restart = new EzVarSequence("Starting point:");
        restart.setNoSequenceSelection();
        channelRestart = new EzVarChannel("Initialization channel :", restart.getVariable(), false);
        System.out.println("initialize()");
        dataSize = new EzVarText("Data size:");
        outputSize = new EzVarText("Output size:");
        paddingSizeX = new EzVarInteger("padding x:",0, Integer.MAX_VALUE,1);
        paddingSizeY = new EzVarInteger("padding y:",0, Integer.MAX_VALUE,1);
        paddingSizeZ = new EzVarInteger("padding z :",0, Integer.MAX_VALUE,1);

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

                dataSize.setVisible(false);
                outputSize.setVisible(false);

                Sequence seq = data.getValue();
                if (seq != null || (seq != null && seq.isEmpty())) {


                    sizeX =  newValue.getSizeX();
                    sizeY = newValue.getSizeY();
                    sizeZ = newValue.getSizeZ();


                    dataSize.setVisible(true);
                    outputSize.setVisible(true);

                    updatePSFSize();
                    updateImageSize();

                    dataShape = new Shape(sizeX, sizeY, sizeZ);

                    if (debug) {
                        System.out.println("Seq changed:" + sizeX + "  "+ Nx);
                        show(  sequenceToArray( seq,0));
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
                System.out.println("weight:" + weightsMethod.getValue()+".");
                System.out.println("weight:" + newValue+".");
                System.out.println("weight:" + weightOptions[3]+".");
                System.out.println("weight:" + weightOptions[3]==newValue);
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
                DoubleArray dataArray;
                // Preparing parameters and testing input

                Sequence dataSeq = data.getValue();
                if(dataSeq!=null){
                    dataArray =    sequenceToArray(dataSeq, channel.getValue()).toDouble();
                    wgtArray = createWeights(dataArray).toDouble();
                    show(wgtArray,"Weight map");
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

        scale = new EzVarDoubleArrayNative("Aspect ratio of a voxel", scaleDef, 2,true);
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

        showFullObject = new EzButton("Show the full (padded) object", new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(debug){
                    System.out.println("showFullObject");
                }
                Sequence fSeq;
                fSeq = new Sequence("Deconvolved image");
                fSeq.copyMetaDataFrom(data.getValue(), false);
                show(solver.getSolution(),fSeq,"Deconvolved "+ data.getValue().getName() + " mu="+solver.getRegularizationLevel() );
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
                    System.out.println("gain"+ gain.getValue());
                    loadParameters(saveFile.getValue());
                    System.out.println("gain"+ gain.getValue());
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



        EzGroup groupStop1 = new EzGroup("Emergency STOP", stopDec);

        addEzComponent(data);
        addEzComponent(channel);
        addEzComponent(psf);
        addEzComponent(channelpsf);
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
        addEzComponent(groupStop1);

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
    protected void updateOutputSize() {

        String text;
        if (Nz==1){
            text= Nx+"x"+Ny;
            scale.setValue(scaleDef[1]);
        }else{
            scale.setValue(scaleDef[2]);
            text= Nx+"x"+Ny+"x"+Nz;
        }
        outputSize.setValue(text);
        if((1.0*Nx*Ny*Nz)>Math.pow(2, 31)){
            throwError("Padded image is too large (>2^31)");
        }
    }




    @Override
    protected void updatePaddedSize() {

        Nx = FFTUtils.bestDimension(sizeX + paddingSizeX.getValue());
        Ny = FFTUtils.bestDimension(sizeY + paddingSizeY.getValue());
        if ((Nz==1)&&(paddingSizeZ.getValue()==0)){
            outputShape = new Shape(Nx, Ny);
        }else{
            Nz= FFTUtils.bestDimension(sizeZ + paddingSizeZ.getValue());
            outputShape = new Shape(Nx, Ny, Nz);
        }
        updateOutputSize();
        if(debug){
            System.out.println(" UpdatePaddedSize" + paddingSizeX.getValue()  + outputShape.toString());
        }
    }

    protected void updatePSFSize() {
        paddingSizeX.setValue(psfSizeX);
        paddingSizeY.setValue(psfSizeY);
        if( (psfSizeZ==1))
            paddingSizeZ.setValue(  0);
        else
            paddingSizeZ.setValue(  psfSizeZ);

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
            System.out.println("Launch it:"+nbIterDeconv.getValue());
        }
        launch();

    }

    protected void launch() {
        solver = new EdgePreservingDeconvolution();

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
            solver.setInitialSolution(objArray);
        }else{
            if(debug){
                System.out.println("restart seq is null:");
            }
        }

        if(singlePrecision.getValue()){
            wgtArray = createWeights(dataArray.toFloat()).toFloat();
            solver.setForceSinglePrecision(true);
        }else{
            wgtArray = createWeights(dataArray.toDouble()).toDouble();
            solver.setForceSinglePrecision(false);
        }
        cursequence = new Sequence("Current Iterate");
        cursequence.copyMetaDataFrom(dataSeq, false);

        solver.setRelativeTolerance(0.0);
        solver.setUseNewCode(false);
        solver.setObjectShape(outputShape);
        solver.setPSF(psfArray);
        solver.setData(dataArray);
        solver.setWeights(wgtArray);
        solver.setEdgeThreshold(epsilon.getValue());
        solver.setRegularizationLevel(mu.getValue());

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
            solver.setScale(scale.getValue());
        }
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
                if (task == OptimTask.FINAL_X) {
                    break;
                }
            }
            if (task == OptimTask.WARNING) {
                break;
            }
            solver.iterate();
        }
        show(ArrayUtils.crop(solver.getBestSolution().asShapedArray(),dataShape),cursequence,"Deconvolved "+ dataSeq.getName() +  " mu="+solver.getRegularizationLevel());


        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                if(true){
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
                        saveSequence(cursequence, outputPath);
                    }


                    if(saveFile.getValue()!=null){
                        saveParamClicked();
                    }
                }
            }
        });

    }





    //If the user call the stop button
    @Override
    public void stopExecution() {
        run = false;

    }
    /**
     *  set default values of the plugin
     */
    @Override
    protected void setDefaultValue() {
        super.setDefaultValue();
        paddingSizeX.setValue(10);
        paddingSizeY.setValue(10);

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