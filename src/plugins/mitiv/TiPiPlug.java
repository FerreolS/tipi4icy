/**
 *
 */
package plugins.mitiv;

import static plugins.mitiv.io.Icy2TiPi.arrayToSequence;

import java.io.File;

import icy.file.Saver;
import icy.gui.frame.progress.FailedAnnounceFrame;
import icy.sequence.Sequence;
import mitiv.array.ShapedArray;
import mitiv.linalg.shaped.ShapedVector;
import plugins.adufour.ezplug.EzPlug;

/**
 * @author ferreol
 *
 */
public abstract class TiPiPlug extends EzPlug  {

    protected  void show(ShapedVector  arr) {
        show(  arr.asShapedArray(),  null,  "" ) ;
    }
    protected  void show(ShapedVector  arr,  String title ) {
        show(  arr.asShapedArray(),  null,  title ) ;
    }
    protected  void show(ShapedArray  arr,  String title ) {
        show(  arr,  null,  title ) ;
    }
    protected void show(ShapedVector  arr, Sequence sequence, String title ) {
        show(  arr.asShapedArray(),  sequence,  title ) ;
    }
    protected  void show(ShapedArray  arr) {
        show(  arr,  null,  "" ) ;
    }
    protected void show(ShapedArray  arr, Sequence sequence, String title ) {
        if (sequence == null )  {

            sequence = new Sequence();
            if (!isHeadLess()){
                addSequence(sequence);
            }
        }
        sequence.beginUpdate();
        sequence =   arrayToSequence(arr, sequence);

        if( sequence.getFirstViewer() == null){
            if (!isHeadLess()){
                addSequence(sequence);
            }
        }

        sequence.endUpdate();
        sequence.setName(title);



    }

    /**
     * print error message
     *
     * @param s
     * the error message
     */
    protected  void throwError(String s){ //FIXME should be in another class
        if(isHeadLess()){
            throw new IllegalArgumentException(s);
        }
        else{
            new FailedAnnounceFrame(s);

        }

    }


    protected void saveSequence(Sequence seq, String path) //FIXME should be elsewhere
    {

        if(path!=null){
            File f = new File(path);
            if (f.isDirectory())
            {
                f = new File(f.getAbsolutePath() + File.separator + seq + ".tif");
            }
            Saver.save(seq, f, false, false);
        }
    }


}
