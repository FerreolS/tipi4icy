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
 * TiPi plugin class using EzPlug Provides methods to interact with icy
 *
 * @author ferreol
 */
public abstract class TiPiPlug extends EzPlug {

    /**
     * Save a sequence in the file pointed by path
     * 
     * @param seq
     *            sequence to save
     * @param path
     *            path of the file
     */
    protected void saveSequence(final Sequence seq, final String path) {
        if (path != null) {
            File f = new File(path);
            if (f.isDirectory()) {
                f = new File(f.getAbsolutePath() + File.separator + seq + ".tif");
            }
            Saver.save(seq, f, false, false);
        }
    }

    /**
     * Show an array in icy
     * 
     * @param arr
     *            ShapedArray to show
     */
    protected void show(final ShapedArray arr) {
        show(arr, null, "");
    }

    /**
     * Show an array in icy by updating the sequence sequence
     * 
     * @param arr
     *            ShapedArray to show
     * @param sequence
     *            sequence to update
     * @param title
     *            title of the plot
     */
    protected void show(final ShapedArray arr, Sequence sequence, final String title) {
        if (sequence == null) {

            sequence = new Sequence();
            if (!isHeadLess()) {
                addSequence(sequence);
            }
        }
        sequence.beginUpdate();
        sequence = arrayToSequence(arr, sequence);

        if (sequence.getFirstViewer() == null) {
            if (!isHeadLess()) {
                addSequence(sequence);
            }
        }

        sequence.endUpdate();
        sequence.setName(title);

    }

    /**
     * Show an array in icy in a new sequence/window
     * 
     * @param arr
     *            ShapedArray to show
     * @param title
     *            title of the plot
     */
    protected void show(final ShapedArray arr, final String title) {
        show(arr, null, title);
    }

    /**
     * Show a vector in icy in a new sequence/window
     * 
     * @param arr
     *            ShapedVector to show
     */
    protected void show(final ShapedVector arr) {
        show(arr.asShapedArray(), null, "");
    }

    /**
     * Show a vector in icy by updating the sequence sequence
     * 
     * @param arr
     *            ShapedArray to show
     * @param sequence
     *            sequence to update
     * @param title
     *            title of the plot
     */
    protected void show(final ShapedVector arr, final Sequence sequence, final String title) {
        show(arr.asShapedArray(), sequence, title);
    }

    /**
     * Show a vector in icy in a new sequence/window
     * 
     * @param arr
     *            ShapedVector to show
     * @param title
     *            title of the plot
     */
    protected void show(final ShapedVector arr, final String title) {
        show(arr.asShapedArray(), null, title);
    }

    /**
     * print error message
     *
     * @param s
     *            the error message
     */
    protected void throwError(final String s) {
        if (isHeadLess()) {
            throw new IllegalArgumentException(s);
        } else {
            new FailedAnnounceFrame(s);

        }

    }

}
