/*
 * Copyright (c) 2017 Ferréol Soulez ferreol.soulez@univ-lyon1.fr
 *
 * This file is part of microTiPi
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

package plugins.ferreol.demics;

/**
 * Define a class for meta-data
 * @author ferreol
 *
 */
public class MicroscopeMetadata {
    //Just a public object with all psf values inside
    /**
     * Lateral size of a pixel in nm
     */
    public double dx     = 64.5;
    /**
     * Lateral size of a pixel in nm
     */
    public double dy     = 64.5;
    /**
     * axial sampling step size
     */
    public double dz      = 160;
    /**
     * number of pixels along x
     */
    public int    nx     = 256;
    /**
     * number of pixels along  y
     */
    public int    ny     = 256;
    /**
     * number of planes
     */
    public int    nz      = 128;
    /**
     * Numerical aperture
     */
    public double na      = 1.4;
    /**
     * Wavelength in nm
     */
    public double lambda  = 542;
    /**
     * Refractive index of the immersion medium
     */
    public double ni      = 1.518;

    @Override
    public String toString(){
        return new String("dx: "+dx+"dy: "+dy+" dz: "+dz+" nx: "+nx+" ny: "+ny+" nz: "+nz+" na "+na+" lambda "+lambda+" ni "+ni);
    }
}
