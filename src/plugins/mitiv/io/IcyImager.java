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

package plugins.mitiv.io;

import static plugins.mitiv.io.Icy2TiPi.arrayToSequence;

import java.io.File;

import icy.file.Saver;
import icy.main.Icy;
import icy.sequence.Sequence;
import mitiv.array.ShapedArray;
import mitiv.linalg.shaped.ShapedVector;
import mitiv.utils.Imager;

/**
 * @author ferreol
 *
 */

public class IcyImager implements Imager{
    protected boolean headless = false;
    protected String title;
    private Sequence sequence;

    /**
     * @param seq
     * @param headless
     */
    public IcyImager(Sequence seq, boolean headless){
        sequence= seq;
        this.headless = headless;
    }
    /**
     * @param seq
     */
    public IcyImager(Sequence seq){
        this.sequence= seq;
    }
    /**
     * @param headless
     */
    public IcyImager(boolean headless){
        this.headless = headless;
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
     * @param headless
     *              true in headless mode
     */
    public static void show(final ShapedArray arr, Sequence sequence, final String title, boolean headless) {
        if (sequence == null )  {
            sequence = new Sequence();
            if (headless){
                Icy.getMainInterface().addSequence(sequence);
            }
        }
        sequence.beginUpdate();
        sequence =   arrayToSequence(arr, sequence);

        if( (sequence.getFirstViewer() == null)&&(!headless)){
            Icy.getMainInterface().addSequence(sequence);
        }

        sequence.endUpdate();
        sequence.setName(title);

    }

    /* (non-Javadoc)
     * @see mitiv.utils.Imager#show(mitiv.array.ShapedArray)
     */
    @Override
    public void show(ShapedArray arr) {
        show(  arr,  sequence,  title ,headless) ;
    }

    /* (non-Javadoc)
     * @see mitiv.utils.Imager#show(mitiv.linalg.shaped.ShapedVector)
     */
    @Override
    public void show(ShapedVector vec) {
        show(  vec.asShapedArray()) ;
    }

    /* (non-Javadoc)
     * @see mitiv.utils.Imager#show(mitiv.array.ShapedArray, java.lang.String)
     */
    @Override
    public void show(ShapedArray arr, String title) {
        this.title = title;
        show(  arr,  sequence,  title ,headless) ;
    }

    /* (non-Javadoc)
     * @see mitiv.utils.Imager#show(mitiv.linalg.shaped.ShapedVector, java.lang.String)
     */
    @Override
    public void show(ShapedVector vec, String title) {
        show(  vec.asShapedArray(),   title ) ;
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
        this.title = title;
        this.sequence = sequence;
        show(arr.asShapedArray());
    }






    /* (non-Javadoc)
     * @see mitiv.utils.Imager#save(mitiv.array.ShapedArray, java.lang.String)
     */
    @Override
    public void save(ShapedVector vec, String path) {
        save(vec.asShapedArray(), path);
    }

    /* (non-Javadoc)
     * @see mitiv.utils.Imager#save(mitiv.linalg.shaped.ShapedVector, java.lang.String)
     */
    @Override
    public void save(ShapedArray arr, String path) {
        save( arr,  sequence,  path);
    }

    /**
     *
     * Save the  image in a file
     * @param arr
     *        array to save
     * @param sequence
     *        sequence to save.
     *        If arr != null it will be copied in sequence
     * @param path
     *        path of the file
     *
     */
    static public void save(ShapedArray arr, Sequence sequence, String path) {
        {
            if (sequence == null )  {
                if(arr == null)
                    return;
                sequence = new Sequence();
            }
            if(arr != null){
                sequence =   arrayToSequence(arr, sequence);
            }

        }
        if (path != null) {
            File f = new File(path);
            if (f.isDirectory()) {
                f = new File(f.getAbsolutePath() + File.separator + sequence + ".tif");
            }
            Saver.save(sequence, f, false, false);
        }

    }

    /**
     *
     * Save the  image in a file
     * @param sequence
     *        sequence to save.
     * @param path
     *        path of the file
     *
     */
    static public void save(Sequence sequence, String path) {
        {
            if (sequence == null )  {
                return;
            }
        }
        if (path != null) {
            File f = new File(path);
            if (f.isDirectory()) {
                f = new File(f.getAbsolutePath() + File.separator + sequence + ".tif");
            }
            Saver.save(sequence, f, false, false);
        }

    }

}
