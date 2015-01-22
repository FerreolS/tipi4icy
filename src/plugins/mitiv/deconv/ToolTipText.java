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

package plugins.mitiv.deconv;

public class ToolTipText {

    static final String deconvolutionSlider = "Update Mu value";
    static final String sequenceImage = "The image on which we will work";
    static final String sequencePSF = "The PSF associated to the image given";
    static final String sequenceWeigth = "The weight map to possibly ignore or minize errors in the image";
    static final String sequenceVariance = "The variance map";
    static final String sequencePixel = "The binary pixel map representing the pixels that we should ignore";
    
    static final String doubleGrtoll = "Relative gradient tolerance for the convergence";
    static final String doubleMu = "Regularization level";
    static final String doubleEpsilon = "Threshold level";
    static final String doubleMaxIter = "Maximum number of iterations, -1 for no limits";
    static final String doublePadding = "Add zero around the image (2D and 3D)";
    static final String doubleGain = "The gain in e-/level";
    static final String doubleNoise = "The readout noise i.e the RMS in e-/pixel";
    static final String doubleBDecTotalIteration = "The maximum number of loop the algorithm is allowed to do, the higher the potentially longer";
    static final String doubleGrtolPhase = "Relative gradient tolerance for the convergence of the phase coeficients";
    static final String doubleGrtolModulus = "Relative gradient tolerance for the convergence of the modulus coeficients";
    static final String doubleGrtolDefocus = "Relative gradient tolerance for the convergence of the defocus coeficients";
    
    static final String doubleDxy = "Lateral pixel size (nm)";
    static final String doubleDz = "Axial pixel size (nm)";
    static final String doubleNxy = "Number of pixels along XY-axis";
    static final String doubleNz = "Number of pixels along Z-axis";
    static final String doubleNa = "Numerical Aperture";
    static final String doubleLambda = "Wavelength (nm)";
    static final String doubleNi = "Refractive index of the immersion medium";
    static final String doubleNalpha = "Number of zernike describing the pupil phase";
    static final String doubleNbeta = "Number of zernike describing the pupil modulus";
    //static final String booleanPSFSplitted = "If the PSF is not centered check this box";
    static final String booleanPSFSplitted = "<html><pre>"
            + "Is the PSF:<br/>"
            + " ---------      ---------     <br/>"
            + " |       |      |0     0|     <br/>"
            + " |   0   |  OR  |       |     <br/>"
            + " |       |      |0     0|     <br/>"
            + " ---------      ---------     <br/>"
            + "</pre></html>";
    static final String booleanRestart = "Restart from previous result, if enabled will start with last image and PSF";
    static final String booleanPositivity = "Limit the negatives values while computing the solution";
    
    static final String textMethod = "Choose the algorithm used to deconvoluate the image";
    static final String textCanal = "Choose the image canal to use for the deconvolution";
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