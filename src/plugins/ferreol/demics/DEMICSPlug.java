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

import icy.gui.frame.progress.AnnounceFrame;
import icy.gui.frame.progress.FailedAnnounceFrame;
import icy.sequence.MetaDataUtil;
import icy.sequence.Sequence;
import icy.util.OMEUtil;
import loci.formats.ome.OMEXMLMetadata;
import loci.formats.ome.OMEXMLMetadataImpl;
import mitiv.array.ShapedArray;
import mitiv.base.Shape;
import mitiv.cost.WeightedData;
import mitiv.utils.FFTUtils;
import mitiv.utils.WeightFactory;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzButton;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarChannel;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarDoubleArrayNative;
import plugins.adufour.ezplug.EzVarFile;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.ezplug.EzVarText;

/**
 * Plugin class for all plugins of DEconvolution MIcroscopy Studio
 * (DEMICS)
 * @author Ferréol
 *
 */
public abstract class DEMICSPlug extends EzPlug  implements Block{

    protected EzVarSequence   data;           // data
    protected EzVarChannel    channel;        // data channel

    protected MicroscopeMetadata meta = null; // metadata of the data


    protected EzVarDouble     logmu, mu;      // deconvolution hyper parameters; mu = 10^(logmu)
    protected EzVarSequence   restart;        // starting point
    protected EzVarChannel    channelRestart; // starting point channel
    protected EzVarBoolean    positivity;     // enforce non negativity
    protected EzButton        startDec, showFullObject;


    protected EzVarText       dataSize;       //
    protected EzVarText       outputSize;     // size of the object after padding
    // optical parameters
    protected EzVarDouble     dy_nm,dx_nm, dz_nm;  //  pixels size in (x,y) and z
    protected EzVarDouble     na  ;             //  numerical aperture
    protected EzVarDouble     lambda;         //  wavelength
    protected EzVarDouble     ni;             //  refractive index of the immersion index

    protected EzVarInteger    nbIterDeconv;   // number of iteration for the deconvolution stage
    protected EzVarBoolean    singlePrecision;// compute in single precision
    protected EzVarDoubleArrayNative scale;   // scale of a voxel should be [1 1 dz/dxy]

    protected EzVarInteger    paddingSizeXY, paddingSizeZ; // number of pixels added in each direction

    // Main variables for the deconvolution part
    protected int sizeX=128, sizeY=128, sizeZ=64; // Input sequence sizes
    protected  int Nx=128,Ny=128, Nz=64;             // Output (padded sequence size)
    protected Shape psfShape = new Shape(Nx, Ny, Nz);
    protected Shape outputShape;
    protected Sequence dataSeq;
    protected Sequence cursequence; // Sequence containing the current solution
    protected Shape dataShape;
    protected ShapedArray wgtArray, dataArray, psfArray, objArray;


    protected EzVarText       weightsMethod;  // Combobox for variance estimation
    protected final String[] weightOptions = new String[]{"None","Inverse covariance map","Variance map","Computed variance"};
    protected EzVarDouble     gain, noise;    // gain of the detector in e-/lvl and detector noise in e-
    protected EzVarSequence weights, deadPixel; // maps of inverse variance and bad pixels
    protected EzButton        showWeight;

    protected EzVarFile       saveFile, loadFile;// xml files to save and load parameters
    protected EzVarBoolean    showIteration;  // show object update at each iteration
    protected EzVarSequence   outputHeadlessImage=null;
    protected EzVarSequence   outputHeadlessWght=null;


    protected String outputPath=null;

    /*********************************/
    /**            DEBUG            **/
    /*********************************/
    private boolean debug = false;      // Show psf steps
    @SuppressWarnings("unused")
    private boolean verbose = false;    // Show some values, need debug to true

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
     * @param badArray
     * @return The weights.
     */

    protected ShapedArray createWeights(ShapedArray datArray, ShapedArray badArray) {
        ShapedArray wgtArray = null;
        Sequence seq;
        WeightedData wd = new WeightedData(datArray);

        if (weightsMethod.getValue() == weightOptions[1]) {
            // A map of weights is provided.
            if ((seq = weights.getValue()) != null) {
                wgtArray =  sequenceToArray(seq);
                if (!wgtArray.equals(datArray)){
                    throwError("Weight map must have the same size than the data");
                }
                wd.setWeights(wgtArray);
            }
        } else if (weightsMethod.getValue() == weightOptions[2]) {
            // A variance map is provided. FIXME: check shape and values.
            if ((seq = weights.getValue()) != null) {
                ShapedArray varArray =  sequenceToArray(seq);
                if (!varArray.equals(datArray)){
                    throwError("Variance map must have the same size than the data");
                }
                wgtArray = WeightFactory.computeWeightsFromVariance(varArray);
                wd.setWeights(wgtArray);
            }
        } else if (weightsMethod.getValue() == weightOptions[3]) {
            // Weights are computed given the gain and the readout noise of the detector.
            double gamma = gain.getValue();
            double sigma = noise.getValue();
            double alpha = 1/gamma;
            double beta = (sigma/gamma)*(sigma/gamma);
            wd.computeWeightsFromData(alpha, beta);
        }
        /*     if ((seq = deadPixel.getValue()) != null) {
            // Account for bad data.
            ShapedArray badArr =  sequenceToArray(seq);
            wd.markBadData(badArr);
        }*/
        if (badArray != null) {
            // Account for bad data.
            if (!badArray.equals(datArray)){
                throwError("Bad data map must have the same size than the data");
            }
            wd.markBadData(badArray);
        }
        return wd.getWeights().asShapedArray();

    }

    /**
     * Function triggered when the data sequence change.
     * Change parameters according to metadata of the data sequence
     *
     */
    protected void dataChanged() {
        dataSize.setVisible(false);
        outputSize.setVisible(false);
        dataSeq = data.getValue();
        if(dataSeq!=null){
            sizeX = dataSeq.getSizeX();
            sizeY = dataSeq.getSizeY();
            sizeZ = dataSeq.getSizeZ();
            // setting restart value to the current sequence
            restart.setValue(dataSeq);
            channelRestart.setValue(channel.getValue());
            updatePaddedSize();
            updateOutputSize();
            updateImageSize();
            dataShape = new Shape(sizeX, sizeY, sizeZ);
            dataSize.setVisible(true);
            outputSize.setVisible(true);

            startDec.setEnabled(true);

            meta = getMetaData(dataSeq);
            dx_nm.setValue(    meta.dx);
            dy_nm.setValue(    meta.dy);
            dz_nm.setValue(     meta.dz);
            if (sizeZ==1) {
                scale.setValue(new double[]{1.0 , dx_nm.getValue()/ dy_nm.getValue()} );
            }else {
                scale.setValue(new double[]{1.0 , dx_nm.getValue()/ dy_nm.getValue(), dx_nm.getValue()/ dz_nm.getValue() } );
            }
        }
    }



    /* (non-Javadoc)
     * The input variable for the protocol
     * @see plugins.adufour.blocks.lang.Block#declareInput(plugins.adufour.blocks.util.VarList)
     */
    @Override
    public void declareInput(VarList inputMap) {
        inputMap.add("image", data.getVariable());
        inputMap.add("image channel", channel.getVariable());
        inputMap.add("starting point", restart.getVariable());
        channelRestart = new EzVarChannel("Initialization channel :", restart.getVariable(), false);

        inputMap.add("starting point channel", channelRestart.getVariable());

        inputMap.add("weights Method",weightsMethod.getVariable());
        inputMap.add("deadPixel", deadPixel.getVariable());
        inputMap.add("gain", gain.getVariable());
        inputMap.add("noise", noise.getVariable());

        inputMap.add("mu", mu.getVariable());
        inputMap.add("scale", scale.getVariable());

        inputMap.add("nbIteration", nbIterDeconv.getVariable());
        inputMap.add("positivity", positivity.getVariable());
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
        outputMap.add("outputSize", outputSize.getVariable());
        outputMap.add("output", outputHeadlessImage.getVariable());
        outputMap.add("weightmap", outputHeadlessWght.getVariable());
    }


    /**
     *  set default values of the plugin
     */
    protected void setDefaultValue() {
        weightsMethod.setValue( weightOptions[3]);
        data.setNoSequenceSelection();
        deadPixel.setNoSequenceSelection();


        paddingSizeZ.setValue(30);
        deadPixel.setNoSequenceSelection();


        if (!isHeadLess()) {
            outputSize.setEnabled(false);
            dataSize.setEnabled(false);
            mu.setEnabled(false);
        }
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
        dataSize.setValue(text);
    }

    /**
     * update metaData of the data sequence using the value indicated in the plugin
     */
    protected void updateMetaData() {
        Sequence seq = data.getValue();
        if (seq != null) {
            ome.xml.meta.OMEXMLMetadata newMetdat = MetaDataUtil.generateMetaData(seq, false);
            newMetdat.setPixelsPhysicalSizeX(OMEUtil.getLength(dx_nm.getValue()*1E-3), 0);
            newMetdat.setPixelsPhysicalSizeY(OMEUtil.getLength(dy_nm.getValue()*1E-3), 0);
            newMetdat.setPixelsPhysicalSizeZ(OMEUtil.getLength(dz_nm.getValue()*1E-3), 0);
            seq.setMetaData((OMEXMLMetadataImpl) newMetdat); //FIXME may not working now
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
        outputSize.setValue(text);
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
        }
        if (metDat.getInstrumentCount() > 0) {
            try {
                meta.na      = metDat.getObjectiveLensNA(0, 0);
                //meta.lambda  = metDat.getChannelEmissionWavelength(0, 0).getValue().doubleValue()*1E6;  //I suppose the value I will get is in um
            } catch(Exception e){
                System.out.println("Failed to get some metadatas, will use default values for na, lambda");
            }
        }

        //If no instrument found, at least we have the right image size
        meta.nz      = seq.getSizeZ();
        meta.dx     = seq.getPixelSizeX()*1E3;
        meta.dy     = seq.getPixelSizeY()*1E3;
        meta.dz      = seq.getPixelSizeZ()*1E3;
        if (na !=null) {
            meta.na      = na.getValue();
        }
        if (lambda !=null) {
            meta.lambda  = lambda.getValue();
        }
        if (ni !=null) {
            meta.ni      = ni.getValue();
        }
        return meta;
    }

    /**
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
}
