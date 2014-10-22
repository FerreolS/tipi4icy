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

import java.awt.image.BufferedImage;
import mitiv.linalg.LinearConjugateGradient;

public class ThreadCG extends Thread {

    MitivDeconvolution deconv;

    boolean stop = false;
    boolean hasjob = false;

    int nextJobValue;
    int nextJobJob;

    private boolean compute3D = false;
    
    public ThreadCG(MitivDeconvolution deconv){
        this.deconv = deconv;
    }
    
    public void compute3D(){
        this.compute3D = true;
    }

    public void prepareNextJob(int tmp, int job){
        this.nextJobValue = tmp;
        this.nextJobJob = job;
        hasjob = true;
    }

    public void run() {
        while (!stop) {
            if (hasjob) {
                hasjob = false; //first because if while computing a new job appear, we will not miss it
                deconv.updateProgressBarMessage("Computing");
                deconv.myseq.beginUpdate();
                if (compute3D) {
                    deconv.nextJob3D(nextJobValue, nextJobJob);
                    deconv.updateProgressBarMessage("Done");
                } else {
                    BufferedImage buffered = deconv.nextJob(nextJobValue, nextJobJob);
                    deconv.updateImage(buffered, nextJobValue);

                    //If we have not finish the computation we will continue it later
                    if (deconv.getOutputValue() == LinearConjugateGradient.IN_PROGRESS) {
                        hasjob = true;
                    }else{
                        deconv.updateProgressBarMessage("Done");
                    }
                }
                deconv.myseq.endUpdate();
            }
            try {
                sleep(50);
            } catch (Exception e) {
                System.err.println("Bad sleep thread");
            }
        }
    }

    public void cancel(){
        stop = true;
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
