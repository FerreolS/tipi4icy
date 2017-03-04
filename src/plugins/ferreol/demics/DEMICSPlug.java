/**
 *
 */
package plugins.ferreol.demics;

import static plugins.mitiv.io.Icy2TiPi.sequenceToArray;

import icy.gui.frame.progress.AnnounceFrame;
import icy.sequence.MetaDataUtil;
import icy.sequence.Sequence;
import icy.util.OMEUtil;
import loci.common.services.ServiceException;
import loci.formats.ome.OMEXMLMetadata;
import loci.formats.ome.OMEXMLMetadataImpl;
import microTiPi.microscopy.MicroscopeMetadata;
import mitiv.array.ShapedArray;
import mitiv.base.Shape;
import mitiv.cost.WeightedData;
import mitiv.utils.FFTUtils;
import mitiv.utils.WeightFactory;
import plugins.adufour.ezplug.EzButton;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarChannel;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarDoubleArrayNative;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.ezplug.EzVarText;
import plugins.mitiv.TiPiPlug;

/**
 * @author ferreol
 *
 */
public abstract class DEMICSPlug extends TiPiPlug {

    protected EzVarSequence   data;           // data
    protected EzVarChannel    channel;        // data channel


    protected EzVarDouble logmu, mu; // deconvolution hyper parameters; mu = 10^(logmu)
    protected EzVarSequence   restart;        // starting point
    protected EzVarChannel    channelRestart; // starting point channel
    protected EzVarBoolean    positivity;     // enforce non negativity
    protected EzButton startDec, stopDec,  showFullObject;


    protected EzVarText       dataSize;       //
    protected EzVarText       outputSize;     // size of the object after padding
    // optical parameters
    protected EzVarDouble     dxy_nm, dz_nm;  //  pixels size in (x,y) and z
    protected EzVarDouble     na;             //  numerical aperture
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
    protected boolean run = true;

    protected MicroscopeMetadata meta = null; // metadata of the data


    protected EzVarText       weightsMethod;  // Combobox for variance estimation
    protected final String[] weightOptions = new String[]{"None","Inverse covariance map","Variance map","Computed variance"};
    protected EzVarDouble     gain, noise;    // gain of the detector in e-/lvl and detector noise in e-
    protected EzVarSequence weights, deadPixel; // maps of inverse variance and bad pixels
    protected EzButton        showWeight;



    /*********************************/
    /**            DEBUG            **/
    /*********************************/
    private boolean debug = false;      // Show psf steps
    private boolean verbose = false;    // Show some values, need debug to true

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
    protected MicroscopeMetadata getMetaData(Sequence seq){ // Should be elsewhere
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

    protected void updateMetaData() {
        Sequence seq = data.getValue();
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



    protected void updateImageSize() {
        String text ;
        if (Nz==1)
            text= sizeX+"x"+sizeY;
        else
            text= sizeX+"x"+sizeY+"x"+sizeZ;
        dataSize.setValue(text);
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

    protected ShapedArray createWeights(ShapedArray datArray) { //FIXME should be elsewhere
        ShapedArray wgtArray = null;
        Sequence seq;
        WeightedData wd = new WeightedData(datArray);

        if (weightsMethod.getValue() == weightOptions[1]) {
            // A map of weights is provided.
            if ((seq = weights.getValue()) != null) {
                wgtArray =  sequenceToArray(seq);
                wd.setWeights(wgtArray);
            }
        } else if (weightsMethod.getValue() == weightOptions[2]) {
            // A variance map is provided. FIXME: check shape and values.
            if ((seq = weights.getValue()) != null) {
                ShapedArray varArray =  sequenceToArray(seq);
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
        if ((seq = deadPixel.getValue()) != null) {
            // Account for bad data.
            ShapedArray badArr =  sequenceToArray(seq);
            wd.markBadData(badArr);
        }
        return wd.getWeights().asShapedArray();

    }


    /**
     *
     */
    protected void dataChanged() {
        dataSeq = data.getValue();
        if(dataSeq!=null){
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

    /**
     * Print the size of the deconvolved image  in the plugin
     */
    protected void updateOutputSize() {
        String text = Nx+"x"+Ny+"x"+Nz;
        outputSize.setValue(text);
        if((1.0*Nx*Ny*Nz)>Math.pow(2, 30)){
            throwError("Padded image is too large (>2^30)");
        }
    }


    /**
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
}
