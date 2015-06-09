package plugins.mitiv.io;

import java.io.File;

import mitiv.array.Double2D;
import mitiv.array.Double3D;
import mitiv.array.ShapedArray;
import mitiv.io.MdaFormat;
import icy.gui.frame.progress.AnnounceFrame;
import icy.sequence.Sequence;
import icy.type.collection.array.Array1DUtil;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzVar;
import plugins.adufour.ezplug.EzVarFile;
import plugins.adufour.ezplug.EzVarListener;
import plugins.adufour.ezplug.EzVarSequence;

/**
 * MiTivGlobalDeconv is a blind deconvolution tool built on the same basis than
 * MiTivTotalVariation. The blind deconvolution process is trying to guess the PSF 
 * and then the it is a standard deconvolution. By iterating on the result we can affine 
 * the PSF until the result is good enough for the user.
 * 
 * @author light
 *
 */
public class ImageToMda extends EzPlug{

    String defaultPath = "";

    EzVarSequence input = new EzVarSequence("Input");
    EzVarFile text = new EzVarFile("Output",defaultPath);


    @Override
    public void clean() {
        // TODO Auto-generated method stub

    }

    @Override
    protected void execute() {
        Sequence in = input.getValue();
        try {
            if (in == null) {
                throw new IllegalArgumentException("No input specified");
            }
            if (in.getSizeT() != 1) {
                throw new IllegalArgumentException("No 4D data input");
            }
            int[] dims;
            if (in.getSizeZ() != 1) {
                dims = new int[]{in.getSizeX(),in.getSizeY(),in.getSizeZ()};
            } else {
                dims = new int[]{in.getSizeX(),in.getSizeY()};
            }
            double[] inputData = Array1DUtil.arrayToDoubleArray(in.getDataCopyXYZT( 0 ), in.isSignedDataType());
            ShapedArray out;
            if (dims.length == 3) {
                out = Double3D.wrap(inputData, dims);
            } else {
                out = Double2D.wrap(inputData, dims);
            }
            MdaFormat.save(out, text.getValue().getAbsolutePath());
        } catch(Exception e){
            new AnnounceFrame("Oops, Error: "+ e);
            e.printStackTrace();
        }

    }

    @Override
    protected void initialize() {
        addEzComponent(input);
        addEzComponent(text);
        input.addVarChangeListener(new EzVarListener<Sequence>() {

            @Override
            public void variableChanged(EzVar<Sequence> source, Sequence newValue) {
                if (newValue != null && text.getValue() != null) {
                    text.setValue(new File(newValue.getFilename()+".mda"));
                }

            }
        });
    }

}
