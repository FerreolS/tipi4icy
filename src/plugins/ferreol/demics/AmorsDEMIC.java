package plugins.ferreol.demics;

/*
 /*  Copyright (C) 2023  Ferreol Soulez ferreol.soulez@univ-lyon1.fr
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

import static plugins.mitiv.io.Icy2TiPi.arrayToSequence;
import static plugins.mitiv.io.Icy2TiPi.sequenceToArray;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.SwingUtilities;

import icy.file.Loader;
import icy.gui.frame.progress.AnnounceFrame;
import icy.main.Icy;
import icy.plugin.PluginLauncher;
import icy.plugin.PluginLoader;
import icy.sequence.Sequence;
import mitiv.array.ArrayUtils;
import mitiv.base.Shape;
import mitiv.base.indexing.BoundaryConditions;
import mitiv.conv.WeightedConvolutionCost;
import mitiv.cost.DifferentiableCostFunction;
import mitiv.cost.FiniteDifferenceOperator;
import mitiv.cost.HomogeneousHyperbolicTotalVariation;
import mitiv.cost.QuadraticCost;
import mitiv.jobs.DeconvolutionJob;
import mitiv.linalg.shaped.DoubleShapedVectorSpace;
import mitiv.linalg.shaped.FloatShapedVectorSpace;
import mitiv.jobs.AmorsJob;
import mitiv.utils.FFTUtils;
import mitiv.utils.HistoMap;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzButton;
import plugins.adufour.ezplug.EzGroup;
import plugins.adufour.ezplug.EzVar;
import plugins.adufour.ezplug.EzVarChannel;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarListener;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.mitiv.io.DeconvHook;
import plugins.mitiv.io.IcyImager;


/**
 * SimpleDEMIC implements a TV regularized multi-dimensional deconvolution.
 * @author Ferréol
 *
 */
public class AmorsDEMIC extends DEMICSPlug {


    /***************************************************/
    /**         Plugin interface variables            **/
    /***************************************************/

    private EzVarSequence psf; // Point Spread Function

    private EzVarInteger    paddingSizeX,paddingSizeY;

    /** deconvolution tab: **/
    private EzVarDouble  epsilon;
    private EzVarChannel channelpsf;
    private EzGroup ezDeconvolutionGroup2;

	private EzVarInteger    totalNbOfBlindDecLoop;// number of outer loop



    private IcyImager curImager;
	private Sequence PSFsequence;

	private IcyImager PSFImager;




    /****************************************************/
    /**                 VARIABLES                      **/
    /****************************************************/

    static boolean debug =true;

    private int psfSizeX=1,psfSizeY=1,psfSizeZ=1;
    protected DeconvolutionJob PSFdeconvolver =null ;
    protected DifferentiableCostFunction PSFprior =null;

	private AmorsJob amors =null;




    /*********************************/
    /**      Initialization         **/
    /*********************************/
    @Override
    protected void initialize() {
        super.initialize();

        psf = new EzVarSequence("PSF:");
        channelpsf = new EzVarChannel("PSF channel :", psf.getVariable(), false);
        psf.setNoSequenceSelection();

        paddingSizeX = new EzVarInteger("padding x:",0, Integer.MAX_VALUE,1);
        paddingSizeY = new EzVarInteger("padding y:",0, Integer.MAX_VALUE,1);

        dataSizeTxt.setVisible(false);
        outputSizeTxt.setVisible(false);

        ezPaddingGroup = new EzGroup("Padding",paddingSizeX,paddingSizeY,paddingSizeZ);
        ezPaddingGroup.setFoldedState(true);

        paddingSizeX.addVarChangeListener(zeroPadActionListener);
        paddingSizeY.addVarChangeListener(zeroPadActionListener);
        paddingSizeZ.addVarChangeListener(zeroPadActionListener);



        dataEV.addVarChangeListener(new EzVarListener<Sequence>() {
            @Override
            public void variableChanged(EzVar<Sequence> source,
                    Sequence newValue) {
                if (debug) {
                    System.out.println("Seq ch..."+dataEV.getValue());
                }
                modelArray = null;
                dataChanged() ;
                updatePSFSize();

            }

        });


        psf.addVarChangeListener(new EzVarListener<Sequence>() {
            @Override
            public void variableChanged(EzVar<Sequence> source,
                    Sequence newValue) {
                if (debug) {
                    System.out.println("PSF changed"+psf.getValue());
                }
                if (newValue != null){
					if(newValue.isEmpty()) {
                    	psfSizeX = Math.max(1,newValue.getSizeX());
                    	psfSizeY =  Math.max(1,newValue.getSizeY());
                    	psfSizeZ =  Math.max(1,newValue.getSizeZ());

	                    modelArray = null;
    	                updatePSFSize();
        	            updateImageSize();


            	        if (debug) {
                	        System.out.println("PSF changed:" + psfSizeX + "  "+ psfSizeY);
                    	}
					}
                }
            }

        });




        /****************************************************/
        /**                    DECONV GROUP                  **/
        /****************************************************/
        epsilon = new EzVarDouble("Threshold level:",1E-2,0.,Double.MAX_VALUE,1.0);

        ezDeconvolutionGroup2 = new EzGroup("More  parameters",epsilon,scale,positivityEV,singlePrecision);
        ezDeconvolutionGroup2.setFoldedState(true);

		totalNbOfBlindDecLoop = new EzVarInteger(  "Number of outer loops:",2,0,Integer.MAX_VALUE ,1);
        nbIterDeconv = new EzVarInteger("Number of inner iterations: ",50,1,Integer.MAX_VALUE ,1);
        ezDeconvolutionGroup = new EzGroup("Deconvolution parameters",totalNbOfBlindDecLoop,nbIterDeconv,logmu,mu,ezDeconvolutionGroup2);

        startDecButton = new EzButton("Start Blind Deconvolution", new ActionListener() {
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




        addEzComponent(dataEV);
        addEzComponent(channelEV);
        addEzComponent(psf);
        addEzComponent(channelpsf);
        addEzComponent(restartEV);
        addEzComponent(channelRestartEV);
        addEzComponent(dataSizeTxt);
        addEzComponent(outputSizeTxt);
        addEzComponent(ezPaddingGroup);

        ezWeightingGroup = new EzGroup("Weighting",weightsMethod,weightsSeq,gain,noise,badpixMap,showWeightButton);
        ezWeightingGroup.setFoldedState(true);
        addEzComponent(ezWeightingGroup);
        addEzComponent(ezDeconvolutionGroup);
        addEzComponent(startDecButton);
        addEzComponent(groupVisu);
        addEzComponent(saveFile);
        addEzComponent(saveParam);
        addEzComponent(loadParam);

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
            outputHeadlessWght = new EzVarSequence("Computed weight");
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

        try {
            startDecButton.setText("Emergency stop");

            dataSeq = dataEV.getValue();
            Sequence psfSeq = psf.getValue();
            Sequence restartSeq = restartEV.getValue();


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
			
			if (singlePrecision.getValue()) {
				dataArray =  sequenceToArray(dataSeq, channelEV.getValue()).toFloat();
				psfArray =  sequenceToArray(psfSeq,  channelpsf.getValue()).toFloat();
				
				if (restartEV.getValue() != null && restartSeq != null){
					objArray =  sequenceToArray(restartSeq, channelRestartEV.getValue()).toFloat();
					if(debug){
						System.out.println("restart seq:" +restartSeq.getName());
					}
				}else{
					objArray = sequenceToArray(dataSeq, channelEV.getValue()).toFloat();
					if(debug){
						System.out.println("restart seq is null:");
					}
				}
			}else {
				dataArray =  sequenceToArray(dataSeq, channelEV.getValue()).toDouble();
				psfArray =  sequenceToArray(psfSeq,  channelpsf.getValue()).toDouble();
				
				if (restartEV.getValue() != null && restartSeq != null){
					objArray =  sequenceToArray(restartSeq, channelRestartEV.getValue()).toDouble();
					if(debug){
						System.out.println("restart seq:" +restartSeq.getName());
					}
				}else{
					objArray = sequenceToArray(dataSeq, channelEV.getValue()).toDouble();
					if(debug){
						System.out.println("restart seq is null:");
					}
				}
			}
            createWeights(true);



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

            cursequence = new Sequence("Current object");
            cursequence.copyMetaDataFrom(dataSeq, false);
            curImager = new IcyImager(cursequence, isHeadLess());

            DeconvHook ObjectHook = new DeconvHook(curImager, dataShape,null, false);
            DeconvHook ObjectHookfinal = new DeconvHook(curImager, dataShape,"Deconvolved "+dataSeq.getName(), false);
            
			
            PSFsequence = new Sequence("Current PSF");
            PSFsequence.copyMetaDataFrom(dataSeq, false);
            PSFImager = new IcyImager(PSFsequence, isHeadLess());

            DeconvHook PSFHook = new DeconvHook(PSFImager, dataShape,null, false);
            DeconvHook PSFHookfinal = new DeconvHook(PSFImager, dataShape,"Estimated PSF "+dataSeq.getName(), false);

            buildVectorSpaces();
			psfArray = ArrayUtils.extract(psfArray, outputShape,0.0); //Padding to the right size


            fdata =  WeightedConvolutionCost.build( objectSpace, dataSpace);
            fdata.setData(dataArray);
            fdata.setWeights(wgtArray,true);
            fdata.setPSF(psfArray);            
			
			objArray = ArrayUtils.extract(objArray, outputShape, fdata.getWeightedMean()); //Padding to the right size



           	//fprior = new HomogeneousHyperbolicTotalVariation(objectSpace, epsilon.getValue(), scale.getValue());
            if  (singlePrecision.getValue()){
                fprior = new QuadraticCost(new FiniteDifferenceOperator((FloatShapedVectorSpace)objectSpace));
            }else{
                fprior = new QuadraticCost(new FiniteDifferenceOperator((DoubleShapedVectorSpace) objectSpace));
            }
            deconvolver  = new DeconvolutionJob( fdata,  mu.getValue(),fprior,  positivityEV.getValue(),nbIterDeconv.getValue(),  ObjectHook,  ObjectHookfinal);
			curImager.show(objArray, "obj");

           	//PSFprior = new HomogeneousHyperbolicTotalVariation(objectSpace, epsilon.getValue(), scale.getValue());
			
            
            if  (singlePrecision.getValue()){
                PSFprior = new QuadraticCost(new FiniteDifferenceOperator((FloatShapedVectorSpace)objectSpace,BoundaryConditions.PERIODIC));
            }else{
                PSFprior = new QuadraticCost(new FiniteDifferenceOperator((DoubleShapedVectorSpace) objectSpace, BoundaryConditions.PERIODIC));
            }
            
            //PSFprior = new QuadraticCost(objectSpace);
			PSFdeconvolver  = new DeconvolutionJob( fdata,  1.0,PSFprior,  true, nbIterDeconv.getValue(),  PSFHook,  PSFHookfinal);
			PSFImager.show(psfArray, "PSF");

			amors = new AmorsJob(totalNbOfBlindDecLoop.getValue(), deconvolver,PSFdeconvolver,null, debug);

            amors.blindDeconv(objArray,psfArray);



            if (weightsMethod.getValue() == weightOptions[4]) {
                modelArray =  fdata.getModel(objArray).asShapedArray();
                HistoMap hm = new HistoMap(modelArray, dataArray, badpixArray);
                gain.setValue(hm.getAlpha());
                noise.setValue(Math.sqrt(hm.getBeta())/hm.getAlpha());
            }


            SwingUtilities.invokeLater(new Runnable() {
                @Override
                public void run() {
                    if(debug){
                        System.out.println("invoke later");
                    }
                    restartEV.setValue(cursequence);
                    channelRestartEV.setValue(0);

                    if (isHeadLess()) {
                        if(outputHeadlessImage==null){
                            outputHeadlessImage  = new EzVarSequence("Output Image");
                        }
                        outputHeadlessImage.setValue(cursequence);
                        if (outputHeadlessWght==null) {
                            outputHeadlessWght = new EzVarSequence("Computed weightsSeq");
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
            startDecButton.setText("Start Deconvolution");
        }
    }



    //If the user call the stop button
    @Override
    public void stopExecution() {
        if (amors !=null)
            amors.abort();
    }
    /**
     *  set default values of the plugin
     */
    @Override
    protected void setDefaultValue() {
        super.setDefaultValue();
		weightsMethod.setValue( weightOptions[0]);
        paddingSizeX.setValue(0);
        paddingSizeY.setValue(0);
        paddingSizeZ.setValue(0);
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
                    dataEV.setValue(Loader.loadSequence(args[i+1], 0, false));
                    dataEV.valueChanged(dataEV.getVariable(), null, dataEV.getValue());
                    if(i+3 >= args.length)
                        break;
                    if( args[i+1].startsWith("-"))
                        break;
                    if(args[i+2].equalsIgnoreCase("-c")){
                        channelEV.setValue(Integer.parseInt(args[i+3]));
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
                    outputPath = args[i+1];
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
                    System.out.println("Usage: ");
                    System.out.println("java -jar -Xms24G icy.jar -hl -x plugins.ferreol.demics.SimpleDEMIC ParametersFile.xml -i DataFile   -c CHANNELNUMBER   -p PSFFile  -c CHANNELNUMBER  -r InitialGuessFile -c CHANNELNUMBER  -o DeconvolvedOutputFile");
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
        if (dataSeq != null) {
            sizeX = dataSeq.getSizeX();
            sizeY = dataSeq.getSizeY();
            sizeZ = dataSeq.getSizeZ();
            updatePaddedSize();
            updateOutputSize();
            updateImageSize();
            dataShape = new Shape(sizeX, sizeY, sizeZ);
        }

    }


    @Override
    protected void saveParamClicked() {
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

    /**
     * Only for test purpose.
     */
    public static void main(final String[] args) {
        // Launch the application.
        Icy.main(args);

        /*
         * Programmatically launch a plugin, as if the user had clicked its
         * button.
         */
        PluginLauncher.start(PluginLoader.getPlugin(AmorsDEMIC.class.getName()));
    }


}
