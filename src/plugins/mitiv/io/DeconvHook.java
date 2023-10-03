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

import mitiv.array.ArrayUtils;
import mitiv.base.Shape;
import mitiv.invpb.Deconvolution;
import mitiv.invpb.SmoothInverseProblem;
import mitiv.utils.Imager;
import mitiv.utils.TiPiHook;

/**
 * Hook to plot current object during  deconvolution
 * @author ferreol
 *
 */
public class DeconvHook implements TiPiHook{

    private Imager curImager;
    private Shape outShape;
    private boolean debug;
    private String title=null;

    /**
     * @param imager    An abstract class for plot/save
     * @param outShape  shape of the shown object (depending on padding)
     * @param title     title of the plot
     * @param debug
     *
     */
    public DeconvHook(Imager imager, Shape outShape, String title, boolean debug) {
        this.curImager = imager;
        this.outShape = outShape;
        this.debug = debug;
        this.title = title;
    }
    /* (non-Javadoc)
     * @see mitiv.utils.TiPiHook#run(java.lang.Object, int)
     */
    @Override
    public void run(Object caller, int iter) {
        String titlet;
        if(title==null)
            titlet = "Current mu="+((SmoothInverseProblem) caller).getRegularizationLevel() +"it:"+((SmoothInverseProblem) caller).getIterations();
        else
            titlet = title;
        curImager.show(ArrayUtils.extract(((Deconvolution) caller).getSolution(),outShape),titlet);

        if (debug){
            System.out.println("Cost "+((SmoothInverseProblem) caller).getCost() );
        }
    }


}
