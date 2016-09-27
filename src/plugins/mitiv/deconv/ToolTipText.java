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

/**
 * In this class is contained all the tooltip texts. It can be usefull for correcting mistakes or
 * new languages.
 * 
 * @author light
 *
 */
public class ToolTipText {

    static final String deconvolutionSlider = "Update Mu value";
    static final String sequenceImage = "Stack to be deconvolved";
    static final String sequencePSF = "The PSF associated to the image given";
    static final String sequenceWeigth = "Map of space varying noise variance or precision";
    static final String sequenceVariance = "The variance map";
    static final String sequencePixel = "The binary map with 0 for bad pixels";
    
    static final String doubleGrtoll = "Relative gradient tolerance for the convergence";
    static final String doubleMu = "Hyper-parameter";
    static final String doubleEpsilon = "TV threshold Epsilon";
    static final String doubleMaxIter = "Maximum number of iterations";
    static final String doublePadding = "Pad with X lines";
    static final String doubleGain = "The gain in e-/level";
    static final String doubleNoise = "The readout noise i.e the RMS in e-/pixel";
    static final String doubleBDecTotalIteration = "The maximum number of loops of the algorithm, the higher the potentially longer";
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
            + " |       |      |*     *|     <br/>"
            + " |   *   |  OR  |       |     <br/>"
            + " |       |      |*     *|     <br/>"
            + " ---------      ---------     <br/>"
            + "</pre></html>";
    static final String booleanRestart = "Restart from previous result, if enabled will start with last image and PSF";
    static final String booleanPositivity = "Enforce the positivity of the solution";
    static final String booleanCrop = "Crop the output with the same field of view as the input.";
    static final String textMethod = "Choose the algorithm used to deconvoluate the image";
    static final String textCanal = "Choose the image canal to use for the deconvolution";
    static final String textOutput = "The output size as the nearest power of 2 size";
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