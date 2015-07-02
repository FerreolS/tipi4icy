/*
 * This file is part of TiPi (a Toolkit for Inverse Problems and Imaging)
 * developed by the MitiV project.
 *
 * Copyright (c) 2014 the MiTiV project, http://mitiv.univ-lyon1.fr/
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

package plugins.mitiv.reconstruction;

import icy.sequence.Sequence;
import mitiv.array.ArrayFactory;
import mitiv.array.DoubleArray;
import mitiv.array.ShapedArray;
import mitiv.base.Shape;
import mitiv.cost.CompositeDifferentiableCostFunction;
import mitiv.cost.HyperbolicTotalVariation;
import mitiv.cost.QuadraticCost;
import mitiv.deconv.WeightedConvolutionOperator;
import mitiv.invpb.ReconstructionJob;
import mitiv.invpb.ReconstructionSynchronizer;
import mitiv.invpb.ReconstructionViewer;
import mitiv.linalg.LinearOperator;
import mitiv.linalg.shaped.DoubleShapedVector;
import mitiv.linalg.shaped.DoubleShapedVectorSpace;
import mitiv.linalg.shaped.ShapedLinearOperator;
import mitiv.optim.ArmijoLineSearch;
import mitiv.optim.BLMVM;
import mitiv.optim.BoundProjector;
import mitiv.optim.LBFGS;
import mitiv.optim.LineSearch;
import mitiv.optim.MoreThuenteLineSearch;
import mitiv.optim.NonLinearConjugateGradient;
import mitiv.optim.OptimTask;
import mitiv.optim.ReverseCommunicationOptimizer;
import mitiv.optim.SimpleBounds;
import mitiv.optim.SimpleLowerBound;
import mitiv.optim.SimpleUpperBound;
import mitiv.optim.VMLMB;
import mitiv.utils.Timer;
import mitiv.utils.reconstruction.ReconstructionThreadToken;

public class TotalVariationJobForIcy extends ReconstructionJobForIcy implements ReconstructionJob {
    /*****************************************************************************************/
    private double mu = 10.0;

    private double epsilon = 1.0;

    private double gatol = 0.0;

    private double grtol = 1e-3;

    private int limitedMemorySize = 5;

    private double lowerBound = Double.NEGATIVE_INFINITY;

    private double upperBound = Double.POSITIVE_INFINITY;

    private int maxiter = 200;

    private Shape resultShape;

    private DoubleArray data = null;
    private DoubleArray psf = null;
    private DoubleArray weight = null;
    private double fcost = 0.0;
    private DoubleShapedVector gcost = null;
    private Timer timer = new Timer();

    private ReverseCommunicationOptimizer minimizer = null;
    private ReconstructionViewer viewer = null;
    private ReconstructionSynchronizer synchronizer = null;
    private double[] synchronizedParameters = {0.0, 0.0};
    private boolean[] change = {false, false};

    public ReconstructionViewer getViewer() {
        return viewer;
    }
    public void setViewer(ReconstructionViewer rv) {
        viewer = rv;
    }
    public ReconstructionSynchronizer getSynchronizer() {
        return synchronizer;
    }
    public void createSynchronizer() {
        if (synchronizer == null) {
            synchronizedParameters[0] = mu;
            synchronizedParameters[1] = epsilon;
            synchronizer = new ReconstructionSynchronizer(synchronizedParameters);
        }
    }
    public void deleteSynchronizer() {
        synchronizer = null;
    }
    public DoubleArray getData() {
        return data;
    }
    public void setData(DoubleArray data) {
        this.data = data;
    }
    public DoubleArray getPsf() {
        return psf;
    }
    public void setPsf(DoubleArray psf) {
        this.psf = psf;
    }
    @Override
    public ShapedArray getResult() {
        /* Nothing else to do because the actual result is in a vector
         * which shares the contents of the ShapedArray.  Otherwise,
         * some kind of synchronization is needed. */
        return result;
    }

    public void setResult(ShapedArray result) {
        this.result = result;
    }
    public void setMaximumIterations(int value) {
        maxiter = value;
    }
    public void setLimitedMemorySize(int value) {
        limitedMemorySize = value;
    }
    public void setRegularizationWeight(double value) {
        mu = value;
    }
    public void setRegularizationThreshold(double value) {
        epsilon = value;
    }
    public void setAbsoluteTolerance(double value) {
        gatol = value;
    }
    public void setRelativeTolerance(double value) {
        grtol = value;
    }
    @Override
    public double getRelativeTolerance() {
        return grtol;
    }
    public void setLowerBound(double value) {
        lowerBound = value;
    }
    public void setUpperBound(double value) {
        upperBound = value;
    }
    public void setWeight(DoubleArray W){
        this.weight = W;
    }
    public void setOutputShape(Shape shape){
        resultShape = shape;
    }

    public void setPositivity(boolean bool){
        lowerBound = bool ? 0 : Double.NEGATIVE_INFINITY;
    }

    /********************************************************************************************************************************/

    /**
     *
     * @param sequence
     * @param token
     */
    public TotalVariationJobForIcy(Sequence sequence, ReconstructionThreadToken token) {
        super(sequence, token);
    }

    public TotalVariationJobForIcy(ReconstructionThreadToken token) {
        this(null, token);
    }

    private static void fatal(String reason) {
        throw new IllegalArgumentException(reason);
    }

    @Override
    public void run() {

        //INIT
        timer.start();

        // Check input data and get dimensions.
        if (data == null) {
            fatal("Input data not specified.");
        }
        Shape dataShape = data.getShape();
        Shape psfShape = psf.getShape();
        int rank = data.getRank();

        // Check the PSF.
        if (psf == null) {
            fatal("PSF not specified.");
        }
        if (psf.getRank() != rank) {
            fatal("PSF must have same rank as data.");
        }
        if (resultShape == null) {
            fatal("An output shape must ge given.");
        }

        if (result != null) {
            /* We try to keep the previous result, at least its dimensions
             * must match. */
            for (int k = 0; k < rank; ++k) {
                if (result.getDimension(k) != data.getDimension(k)) {
                    result = null;
                    break;
                }
            }
        }

        // Check the shape of the result.
        for (int k = 0; k < rank; ++k) {
            if (resultShape.dimension(k) < dataShape.dimension(k)) {
                fatal("The dimensions of the result must be at least those of the data.");
            }
            if (resultShape.dimension(k) < psfShape.dimension(k)) {
                fatal("The dimensions of the result must be at least those of the PSF.");
            }

        }

        // Initialize an input and output vector spaces and populate them with
        // workspace vectors.

        DoubleShapedVectorSpace dataSpace = new DoubleShapedVectorSpace(dataShape);
        DoubleShapedVectorSpace resultSpace = new DoubleShapedVectorSpace(resultShape);
        LinearOperator W = null;
        DoubleShapedVector y = dataSpace.create(data);
        DoubleShapedVector x = null;
        if (result != null) {
            x = resultSpace.create(result);
        } else {
            x = resultSpace.create(0.0);
        }
        result = ArrayFactory.wrap(x.getData(), resultShape);

        // Build convolution operator.
        ShapedLinearOperator H = null;
        // FIXME: add a method for that
        WeightedConvolutionOperator A = WeightedConvolutionOperator.build(resultSpace, dataSpace);
        A.setPSF(psf);
        A.setWeights(weight);
        H = A;
        // Build the cost functions
        QuadraticCost fdata = new QuadraticCost(H, y, W);
        HyperbolicTotalVariation fprior = new HyperbolicTotalVariation(resultSpace, epsilon);

        CompositeDifferentiableCostFunction cost = new CompositeDifferentiableCostFunction(1.0, fdata, mu, fprior);
        fcost = 0.0;
        gcost = resultSpace.create();
        timer.stop();

        timer.reset();

        // Initialize the non linear conjugate gradient
        timer.start();
        LineSearch lineSearch = null;
        LBFGS lbfgs = null;
        VMLMB vmlmb = null;
        BLMVM blmvm = null;
        NonLinearConjugateGradient nlcg = null;
        BoundProjector projector = null;
        int bounded = 0;
        if (lowerBound != Double.NEGATIVE_INFINITY) {
            bounded |= 1;
        }
        if (upperBound != Double.POSITIVE_INFINITY) {
            bounded |= 2;
        }
        if (bounded == 0) {
            /* No bounds have been specified. */
            lineSearch = new MoreThuenteLineSearch(0.05, 0.1, 1E-17);
            if (limitedMemorySize > 0) {
                lbfgs = new LBFGS(resultSpace, limitedMemorySize, lineSearch);
                lbfgs.setAbsoluteTolerance(gatol);
                lbfgs.setRelativeTolerance(grtol);
                minimizer = lbfgs;
            } else {
                int method = NonLinearConjugateGradient.DEFAULT_METHOD;
                nlcg = new NonLinearConjugateGradient(resultSpace, method, lineSearch);
                nlcg.setAbsoluteTolerance(gatol);
                nlcg.setRelativeTolerance(grtol);
                minimizer = nlcg;
            }
        } else {
            /* Some bounds have been specified. */
            if (bounded == 1) {
                /* Only a lower bound has been specified. */
                projector = new SimpleLowerBound(resultSpace, lowerBound);
            } else if (bounded == 2) {
                /* Only an upper bound has been specified. */
                projector = new SimpleUpperBound(resultSpace, upperBound);
            } else {
                /* Both a lower and an upper bounds have been specified. */
                projector = new SimpleBounds(resultSpace, lowerBound, upperBound);
            }
            int m = (limitedMemorySize > 1 ? limitedMemorySize : 5);
            //vmlmb = new VMLMB(resultSpace, projector, m, lineSearch);
            //vmlmb.setAbsoluteTolerance(gatol);
            //vmlmb.setRelativeTolerance(grtol);
            //minimizer = vmlmb;
            blmvm = new BLMVM(resultSpace, projector, m);
            blmvm.setAbsoluteTolerance(gatol);
            blmvm.setRelativeTolerance(grtol);
            minimizer = blmvm;
            projector.projectVariables(x);
        }
        timer.stop();
        timer.reset();
        // Launch the non linear conjugate gradient
        OptimTask task = minimizer.start();
        while (!token.isStopped()) {
            if (task == OptimTask.COMPUTE_FG) {
                timer.resume();
                fcost = cost.computeCostAndGradient(1.0, x, gcost, true);
                timer.stop();
            } else if (task == OptimTask.NEW_X || task == OptimTask.FINAL_X) {
                if (viewer != null) {
                    // FIXME: must get values back from the result vector
                    viewer.display(this);
                }
                boolean stop = (task == OptimTask.FINAL_X);
                if (! stop && maxiter >= 0 && minimizer.getIterations() >= maxiter) {
                    System.err.format("Warning: too many iterations (%d).\n", maxiter);
                    stop = true;
                }
                if (stop) {
                    break;
                }
            } else {
                //To see the reason: got to TiPi abstract class ReverseCommunicationOptimizer
                System.err.println("TotalVariationJobForIcy error/warning: "+task+" reason: "+minimizer.getReason());
                break;
            }
            if (synchronizer != null) {
                if (synchronizer.getTask() == ReconstructionSynchronizer.STOP) {
                    break;
                }
                synchronizedParameters[0] = mu;
                synchronizedParameters[1] = epsilon;
                if (synchronizer.updateParameters(synchronizedParameters, change)) {
                    if (change[0]) {
                        mu = synchronizedParameters[0];
                    }
                    if (change[1]) {
                        epsilon = synchronizedParameters[1];
                    }
                }
            }
            task = minimizer.iterate(x, fcost, gcost);
        }
    }
    @Override
    public int getIterations() {
        return (minimizer == null ? 0 : minimizer.getIterations());
    }

    @Override
    public int getEvaluations() {
        return (minimizer == null ? 0 : minimizer.getEvaluations());
    }

    @Override
    public double getCost() {
        return fcost;
    }

    @Override
    public double getGradientNorm2() {
        return (gcost == null ? 0.0 : gcost.norm2());
    }

    public double getLowerBound() {
        return lowerBound;
    }

    @Override
    public double getGradientNorm1() {
        return (gcost == null ? 0.0 : gcost.norm1());
    }

    @Override
    public double getGradientNormInf() {
        return (gcost == null ? 0.0 : gcost.normInf());
    }
}
/*
 * Local Variables:
 * mode: Java
 * tab-width: 8
 * indent-tabs-mode: nil
 * c-basic-offset: 4
 * fill-column: 78
 * coding: utf-8
 * ispell-local-dictionary: "american"
 * End:
 */
