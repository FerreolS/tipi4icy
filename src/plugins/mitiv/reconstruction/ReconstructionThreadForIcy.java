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

public class ReconstructionThreadForIcy extends Thread {

    private ReconstructionThreadToken token;
    private ReconstructionJobForIcy job;

    public ReconstructionThreadForIcy(ReconstructionThreadToken token) {
        this.token = token;
    }

    public void setToken(ReconstructionThreadToken token){
        this.token = token;
    }

    public void setJob(ReconstructionJobForIcy job){
        this.job = job;
    }

    public void run(){
        while (!token.isExiting()) {        //While we should not quit
            while (!token.isRunning()) {    //We wait for a job and for the order to run
                token.waitForStart();
            }
            if (job != null &&!token.isExiting()) { //We run if we don't have to quit
                job.run();
            } else {
                if (!token.isExiting()) {
                    System.err.println("Running command received but no job to run");
                }
                token.jobFinished();
            }
        }
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
