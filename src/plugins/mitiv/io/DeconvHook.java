/**
 *
 */
package plugins.mitiv.io;

import mitiv.array.ArrayUtils;
import mitiv.base.Shape;
import mitiv.jobs.DeconvolutionJob;
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
            titlet = "Current mu="+((DeconvolutionJob) caller).solver.getRegularizationLevel() +"it:"+((DeconvolutionJob) caller).solver.getIterations();
        else
            titlet = title;
        curImager.show(ArrayUtils.crop(((DeconvolutionJob) caller).solver.getSolution(),outShape),titlet);

        if (debug){
            System.out.println("Cost "+((DeconvolutionJob) caller).solver.getCost() );
        }
    }


}
