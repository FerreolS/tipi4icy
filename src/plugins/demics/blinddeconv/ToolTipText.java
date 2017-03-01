/*
 * This file is part of TiPi (a Toolkit for Inverse Problems and Imaging)
 * developed by the MitiV project.
 *
 * Copyright (c) 2014-2016 the MiTiV project, http://mitiv.univ-lyon1.fr/
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

package plugins.demics.blinddeconv;

/**
 * In this class is contained all the tooltip texts. It can be usefull for correcting mistakes or
 * new languages.
 *
 * @author light
 *
 */
public class ToolTipText {

    public static final String deconvolutionSlider = "Update Mu value";
    public static final String sequenceImage = "Stack to be deconvolved";
    public static final String sequencePSF = "The PSF associated to the image given";
    public static final String sequenceWeigth = "Map of space varying noise variance or precision";
    public static final String sequenceVariance = "The variance map";
    public static final String sequencePixel = "The binary map with 0 for bad pixels";

    public static final String doubleGrtoll = "Relative gradient tolerance for the convergence";
    public static final String doubleMu = "Hyper-parameter";
    public static final String doubleEpsilon = "TV threshold Epsilon";
    public static final String doubleMaxIter = "Maximum number of iterations";
    public static final String doublePadding = "Pad with X lines";
    public static final String doubleGain = "The detector gain in electrons per analog digital unit (ADU)";
    public static final String doubleNoise = "The standard deviation of the readout noise in e-/pixel";
    public static final String doubleBDecTotalIteration = "The maximum number of loops of the algorithm, the higher the potentially longer";
    public static final String doubleGrtolPhase = "Relative gradient tolerance for the convergence of the phase coeficients";
    public static final String doubleGrtolModulus = "Relative gradient tolerance for the convergence of the modulus coeficients";
    public static final String doubleGrtolDefocus = "Relative gradient tolerance for the convergence of the defocus coeficients";

    public static final String doubleDxy = "Lateral pixel size (nm)";
    public static final String doubleDz = "Axial pixel size (nm)";
    public static final String doubleNxy = "Number of pixels along XY-axis";
    public static final String doubleNz = "Number of pixels along Z-axis";
    public static final String doubleNa = "Numerical Aperture";
    public static final String doubleLambda = "Wavelength (nm)";
    public static final String doubleNi = "Refractive index of the immersion medium";
    public static final String doubleNalpha = "Number of zernike describing the pupil phase";
    public static final String doubleNbeta = "Number of zernike describing the pupil modulus";
    //static final String booleanPSFSplitted = "If the PSF is not centered check this box";
    public static final String booleanPSFSplitted = "<html><pre>"
            + "Is the PSF:<br/>"
            + " ---------      ---------     <br/>"
            + " |       |      |       |     <br/>"
            + " |   *   |  or  |       |     <br/>"
            + " |       |      |*      |     <br/>"
            + " ---------      ---------     <br/>"
            + "</pre></html>";
    public static final String booleanRestart = "Restart from previous result, if enabled will start with last image and PSF";
    public  static final String booleanPositivity = "Enforce the positivity of the solution";
    public static final String booleanCrop = "Crop the output with the same field of view as the input.";
    public static final String textMethod = "Choose the algorithm used to deconvoluate the image";
    public static final String textCanal = "Choose the image canal to use for the deconvolution";
    public static final String textOutput = "The output size as the nearest power of 2 size";
}
