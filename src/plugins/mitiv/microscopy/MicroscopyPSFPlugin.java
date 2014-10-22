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

package plugins.mitiv.microscopy;

import icy.sequence.Sequence;
import icy.gui.frame.GenericFrame;
import icy.image.IcyBufferedImage;

import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JScrollPane;
import javax.swing.JTextPane;

import mitiv.array.ArrayUtils;
import mitiv.microscopy.*;
import mitiv.utils.MathUtils;
import plugins.adufour.ezplug.EzButton;
import plugins.adufour.ezplug.EzGroup;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarDoubleArrayNative;
import plugins.adufour.ezplug.EzVarInteger;


/**
 * Microscopy PSF
 * 
 * @author Oui oui
 *
 */
public class MicroscopyPSFPlugin extends EzPlug implements EzStoppable
{
    EzVarInteger Nx;
    EzVarInteger Ny;
    EzVarInteger Nz;
    EzVarInteger NZernike;
    EzVarInteger use_depth_scaling;
    EzVarDouble NA;
    EzVarDouble lambda;
    EzVarDouble ni;
    EzVarDouble ns;
    EzVarDouble nzdepth;
    EzVarDouble dxy;
    EzVarDouble dz;
    EzVarDouble deltaX;
    EzVarDouble deltaY;
    EzVarDoubleArrayNative alpha;
    EzVarDoubleArrayNative beta;
    EzVarDouble zdepth;
    EzVarBoolean rho;
    EzVarBoolean phi;
    EzVarBoolean psi;
    
    boolean stopFlag;
    @Override
    protected void initialize()
    {
        /* Initialize variables */
        NA = new EzVarDouble("NA : Numerical Aperture");
        lambda = new EzVarDouble("\u03BB : Wavelength (nm)");
        ni = new EzVarDouble("ni : Refractive index of the immersion medium");
        ns = new EzVarDouble("ns : Refractive index of the specimen");
        dxy = new EzVarDouble("dxy : Lateral pixel size (nm)");
        dz = new EzVarDouble("dz : Axial pixel size (um)");
        Nx = new EzVarInteger("Nx : Number of samples along lateral X-dimension");
        Ny = new EzVarInteger("Ny : Number of samples along lateral Y-dimension");
        Nz = new EzVarInteger("Nz : Number of samples along lateral Z-dimension");
        NZernike = new EzVarInteger("NZernike : Number of zernike modes");
        deltaX = new EzVarDouble("\u03B4x : Defocus along lateral X-dimension");
        deltaY = new EzVarDouble("\u03B4y : Defocus along lateral Y-dimension");
        alpha = new EzVarDoubleArrayNative("\u03B1 : Zernike coefficients of the modulus",
                new double[][] { new double[] { 2.2, 4.4, 6.6 }, new double[] { 1.1, 3.3, 5.5 } }, false);
        beta = new EzVarDoubleArrayNative("\u03B2 : Zernike coefficients of the phase",
                new double[][] { new double[] { 2.2, 4.4, 6.6 }, new double[] { 1.1, 3.3, 5.5 } }, false);
        zdepth = new EzVarDouble("zd : Depth of a light located in the specimen layer (\u03BCm)");
        use_depth_scaling = new EzVarInteger("PSF are centered on the plan with maximum strehl");
        rho = new EzVarBoolean("Show modulus of the pupil", false);
        phi = new EzVarBoolean("Show phi", false);
        psi = new EzVarBoolean("Show psi", false);

        /* Set the default values */
        NA.setValue(1.4);
        lambda.setValue(542.0);
        ni.setValue(1.518);
        ns.setValue(0.0);
        dxy.setValue(64.5);
        dz.setValue(0.16);
        Nx.setValue(256);
        Ny.setValue(256);
        Nz.setValue(64);
        NZernike.setValue(10);
        deltaX.setValue(0.0);
        deltaY.setValue(0.0);
        alpha.setValue(new double[] {0.});
        beta.setValue(new double[] {1.});
        zdepth.setValue(0.);
        use_depth_scaling.setValue(0);

        /* Add to the interface */
        EzGroup parameterGroup = new EzGroup("Enter Widefield Microscope settings", 
                NA, lambda, ni, Nx, Ny, Nz, dxy, dz, ns, NZernike, alpha, beta,
                deltaX, deltaY, zdepth, use_depth_scaling);
        addEzComponent(parameterGroup);
        
        addEzComponent(rho);
        addEzComponent(phi);
        /* Details button */
        EzButton detailsButton = new EzButton("Help", new ActionListener()
        {
            @Override
            public void actionPerformed(ActionEvent e)
            {
                onDetailsClicked();
            }
        });
        
        addEzComponent(detailsButton);
        
        setTimeDisplay(true);
    }

    @Override
    protected void execute()
    {
        MicroscopyModelPSF1D pupil = new MicroscopyModelPSF1D(NA.getValue(),
                lambda.getValue()*1e-9, ni.getValue(), ns.getValue(),
                zdepth.getValue()*1e-6, dxy.getValue()*1e-9, dz.getValue()*1e-6, Nx.getValue(), Ny.getValue(),
                Nz.getValue(), NZernike.getValue(), use_depth_scaling.getValue());

        pupil.computePSF(alpha.getValue(), beta.getValue(),
               deltaX.getValue(), deltaY.getValue(), zdepth.getValue()*1e-6);
        /*
        pupil.setDefocus(new double[] = {};);
        pupil.setPhi(alpha.getValue());
        pupil.setRho(beta.getValue());
        pupil.computePSF();
        */
        double[] PSF_shift = MathUtils.fftShift3D(pupil.getPSF(), Nx.getValue(), Ny.getValue(), Nz.getValue());        
        Sequence psf3DSequence = new Sequence();
        psf3DSequence.setName("PSF");
        for (int k = 0; k < Nz.getValue(); k++)
        {
            psf3DSequence.addImage(new IcyBufferedImage(Nx.getValue(), Ny.getValue(), MathUtils.getArray(PSF_shift, Nx.getValue(), Ny.getValue(), k)));
        }
        MathUtils.stat(pupil.getPSF());

        addSequence(psf3DSequence);
        
        if( rho.getValue() == true )
        {
           Sequence pupilModulusSequence = new Sequence();
           pupilModulusSequence.setName("Modulus");
           //pupilModulusSequence.addImage(ArrayUtils.doubleAsBuffered(MathUtils.fftShift1D(pupil.getRho(), Nx.getValue(), Ny.getValue(), 1),
            //                1, Nx.getValue(), Ny.getValue()));
           pupilModulusSequence.addImage(new IcyBufferedImage(Nx.getValue(), Ny.getValue(), pupil.getRho()));
           addSequence(pupilModulusSequence);
        }
        
        if( phi.getValue() == true )
        {
           Sequence pupilModulusSequence = new Sequence();
           pupilModulusSequence.setName("Phase");
           pupilModulusSequence.addImage(ArrayUtils.doubleAsBuffered(MathUtils.fftShift1D(pupil.getPhi(), Nx.getValue(), Ny.getValue(), 1),
                            1, Nx.getValue(), Ny.getValue()));
           addSequence(pupilModulusSequence);
        }

    }

    @Override
    public void clean() {
        // TODO Auto-generated by Icy4Eclipse
    }
    
    @Override
    public void stopExecution()
    {
        // this method is from the EzStoppable interface
        // if this interface is implemented, a "stop" button is displayed
        // and this method is called when the user hits the "stop" button
        stopFlag = true;
    }
    
    /**
     * Action performed when the user clicks on the details button
     */
    private void onDetailsClicked()
    {
        String    title   = "Noise models";
        JTextPane message = new JTextPane();
        message.setEditable(false);
        message.setContentType("text/html");
        message.setText(
            "<p>" +
                "In what follows, <tt>A</tt> represents the input clean sequence, " +
                "while <tt>B</tt> is the output noisy sequence." +
            "</p>" +
            "<h2>White additive Gaussian noise</h2>" +
            "<p>" +
                "1 parameter: <tt>sigma &gt;= 0</tt> (standard deviation of the Gaussian random variables)." +
            "</p>" +
            "<p>" +
                "Exemple deltaX = ni/lambda/4 = 700185 and deltaY = ni/lambda = 2800740" +
                "zdepth = nz*dz = 10.24 um" +
                "This noise model enforces <tt>B = A + n</tt>, where <tt>n</tt> is a " +
                "random sequence such that the samples <tt>n(x,y,z,t,c)</tt> are random " +
                "independant variables following a Gaussian probability distribution " +
                "of mean <tt>0</tt> and variance <tt>sigma^2</tt>." +
            "</p>" +
            "<h2>Poisson noise</h2>" +
            "<p>" +
                "No parameter." +
            "</p>" +
            "<p>" +
                "In this model, each output sample <tt>B(x,y,z,t,c)</tt> is generated from a Poisson random " +
                "distribution of intensity <tt>A(x,y,z,t,c)</tt>, which is supposed to be <tt>&gt;=0</tt>. " +
                "If <tt>A(x,y,z,t,c) &lt; 0</tt>, then <tt>B(x,y,z,t,c)</tt> is set to <tt>NaN</tt>." +
            "</p>" +
            "<h2>Hints</h2>" +
            "<p>" +
                "<ul>" +
                    "<li>Use the Icy colormap in the lookup table at the left to change the colormap model</li>" +
                    "<li>Use 3D OrthoViewer to see YX and XZ views. To activate this view mode," +
                    "select the Orthoviewer logo in the viewer window (in the dropdown list that defaults to 2D mode).</li>" +
                "</ul>" +
            "</p>"
        );
        Dimension dim = message.getPreferredSize();
        dim.setSize(600, dim.getHeight()+100);
        message.setPreferredSize(dim);
        JScrollPane scroll = new JScrollPane(message);
        dim = scroll.getPreferredSize();
        dim.setSize(600, 500);
        scroll.setPreferredSize(dim);
        GenericFrame infoFrame = new GenericFrame(title, scroll);
        infoFrame.addToMainDesktopPane();
        infoFrame.setVisible(true);
        infoFrame.requestFocus();
    }
    
    /**
     * 
     */
    public double getNa() {
        return NA.getValue();
    }

    /**
     * 
     */
    public double getLambda() {
        return lambda.getValue()*1e-9;
    }

    /**
     * 
     */
    public double getNi() {
        return ni.getValue();
    }

    /**
     
     * 
     */
    public double getDxy() {
        return dxy.getValue()*1e-9;
    }

    /**
     * 
     */
    public double getDz() {
        return dz.getValue()*1e-6;
    }

    /**
     * 
     */
    public double getNx() {
        return Nx.getValue();
    }

    /**
     * 
     */
    public double getNy() {
        return Ny.getValue();
    }

    /**
     * 
     */
    public double getNz() {
        return Nz.getValue();
    }

    /**
     * 
     */
    public double getJ() {
        return NZernike.getValue();
    }

    /**
     * 
     */
    public double getDeltax() {
        return deltaX.getValue();
    }

    /**
     * 
     */
    public double getDeltay() {
        return deltaY.getValue();
    }

    /**
     * 
     */
    public double[] getAlpha() {
        return alpha.getValue();
    }

    /**
     * 
     */
    public double[] getBeta() {
        return beta.getValue();
    }

    /**
     * 
     */
    public boolean getRho() {
        return rho.getValue();
    }

    /**
     * 
     */
    public boolean getPhi() {
        return phi.getValue();
    }

    /**
     * 
     */
    public boolean getPsi() {
        return psi.getValue();
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