package plugins.mitiv.mda;

import java.io.File;
import java.util.ArrayList;

import mitiv.array.ShapedArray;
import mitiv.io.MdaFormat;
import icy.gui.frame.progress.AnnounceFrame;
import icy.image.IcyBufferedImage;
import icy.sequence.Sequence;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzVarFileArray;
import plugins.mitiv.old.io.IcyBufferedImageUtils;

/**
 * MiTivGlobalDeconv is a blind deconvolution tool built on the same basis than
 * MiTivTotalVariation. The blind deconvolution process is trying to guess the PSF 
 * and then the it is a standard deconvolution. By iterating on the result we can affine 
 * the PSF until the result is good enough for the user.
 * 
 * @author light
 *
 */
public class MdaToImage extends EzPlug{

    String defaultPath = "/path/to/input";  
    EzVarFileArray text = new EzVarFileArray("Input",defaultPath);

    @Override
    public void clean() {
        // TODO Auto-generated method stub
    }

    @Override
    protected void execute() {
        File[] files = text.getValue();
        for (int i = 0; i < files.length; i++) {

            try {
                ShapedArray in = MdaFormat.load(files[i].getAbsolutePath());
                ArrayList<IcyBufferedImage> tmp = IcyBufferedImageUtils.arrayToImage(in);
                Sequence out = new Sequence(files[i].getName());
                for (int j = 0; j < tmp.size(); j++) {
                    out.setImage(0, j, tmp.get(j));
                }
                addSequence(out);

            } catch(Exception e){
                new AnnounceFrame("Can not open : "+ files[i].getAbsolutePath());
                e.printStackTrace();
            }
        }
    }

    @Override
    protected void initialize() {
        addEzComponent(text);
    }

}
