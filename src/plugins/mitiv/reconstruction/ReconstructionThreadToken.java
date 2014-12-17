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

import java.util.concurrent.locks.Condition;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

public class ReconstructionThreadToken {
    private double[] values = null;
    private int number = 0;

    private final Lock lock;
    private Condition canRun;
    private Condition jobFinished;

    private boolean run = false;
    private boolean stop = false;
    private boolean exit = false;

    public synchronized void setValue(int i, double value) {
        values[i] = value;
    }

    public ReconstructionThreadToken(double[] initialValues) {
        lock = new ReentrantLock();
        canRun = lock.newCondition();
        jobFinished = lock.newCondition();
        number = (initialValues != null ? initialValues.length : 0);
        values = new double[number];
        for (int i = 0; i < number; ++i) {
            values[i] = initialValues[i];
        }
    }

    public synchronized boolean fetchValues(double[] values, boolean[] changed) {
        boolean anyChange = false;
        for (int i = 0; i < number; ++i) {
            if (values[i] != this.values[i]) {
                anyChange = true;
                changed[i] = true;
                values[i] = this.values[i];
            } else {
                changed[i] = false;
            }
        }
        return anyChange;
    }
    
    public void start(){
        start(true);
    }

    public void start(boolean waitForJob){
        run = true;
        stop = false;
        try {
            lock.lock();
            canRun.signal();
        } finally {
            lock.unlock();
        }
        if (waitForJob) {
            waitForJob();
        }
    }

    public void waitForStart(){
        try {
            lock.lock();
            canRun.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        } finally {
            lock.unlock();
        }
    }

    private void waitForJob(){
        try {
            lock.lock();
            jobFinished.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        } finally {
            lock.unlock();
        }
    }

    public void stop(){
        stop = true;
    }

    public void exit(){
        exit = true;
        start();    //We are unlocking the potentials threads that are waitings to run a job

    }

    public void jobFinished(){
        run = false;
        try {
            lock.lock();
            jobFinished.signal();
        } finally {
            lock.unlock();
        }
    }

    public boolean isRunning() {
        return run;
    }

    public boolean isStopped() {
        return stop;
    }

    public boolean isExiting() {
        return exit;
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
