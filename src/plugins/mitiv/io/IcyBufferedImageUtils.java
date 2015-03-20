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
package plugins.mitiv.io;


import icy.image.IcyBufferedImage;
import icy.sequence.Sequence;
import icy.type.collection.array.Array1DUtil;

import java.awt.image.BufferedImage;
import java.util.ArrayList;

import mitiv.array.Double2D;
import mitiv.array.Double3D;
import mitiv.array.Float2D;
import mitiv.array.Float3D;
import mitiv.array.ShapedArray;
import mitiv.base.Shape;
import mitiv.base.Traits;
import mitiv.io.BufferedImageUtils;
import mitiv.utils.CommonUtils;

public class IcyBufferedImageUtils {

    public static ShapedArray imageToArray(ArrayList<IcyBufferedImage> listImage) {
        int width = listImage.get(0).getWidth();
        int height = listImage.get(0).getHeight();
        int sizeZ = listImage.size();
        int[] shape = new int[]{width, height, sizeZ};
        return imageToArray(listImage,shape[0], shape[1], shape[2]);
    }

    public static ShapedArray imageToArray(ArrayList<IcyBufferedImage> listImage, int[] shape) {
        if (shape.length != 3) {
            throw new IllegalArgumentException("Shape should be of size 3, because input is three dimensionnal data");
        }
        return imageToArray(listImage,shape[0], shape[1], shape[2]);
    }

    public static ShapedArray imageToArray(ArrayList<IcyBufferedImage> listImage, int width, int height, int sizeZ) {
        //First we try if we can convert the data directly by using Icy Capacities
        int[] shape = new int[]{width, height, sizeZ};
        try {
            double[] out;
            out = new double[sizeZ*width*height];
            for (int j = 0; j < sizeZ; j++) {
                double[] tmp = listImage.get(j).getDataCopyCXYAsDouble();
                for (int i = 0; i < tmp.length; i++) {
                    out[i+j*tmp.length] = tmp[i];
                }
            }
            return Double3D.wrap(out, shape);
        } catch (Exception e) {
            //System.err.println("Could not directly convert ICY sequence");
            ArrayList<BufferedImage> list = new ArrayList<BufferedImage>(listImage);
            return BufferedImageUtils.imageToArray(list);
        }
    }

    public static ShapedArray imageToArray(IcyBufferedImage image) {
        int[] shape = new int[]{image.getWidth(), image.getHeight()};
        return imageToArray(image, shape[0], shape[1]);
    }

    public static ShapedArray imageToArray(IcyBufferedImage image, int[] shape) {
        if (shape.length != 2) {
            throw new IllegalArgumentException("Shape should be of size 2, because input is two dimensionnal data");
        }
        return imageToArray(image, shape[0], shape[1]);
    }

    public static ShapedArray imageToArray(IcyBufferedImage image, int width,int height) {
        //First we try if we can convert the data directly by using Icy Capacities
        int[] shape = new int[]{width, height};
        try {
            return Double2D.wrap(image.getDataCopyCXYAsDouble(), shape);
        } catch (Exception e) {
            return BufferedImageUtils.imageToArray(image);
        }
    }

    public static ShapedArray imageToArray(Sequence seq, int canal) {
        if (seq.getSizeZ() == 1) {
            return imageToArray(seq, Shape.make(seq.getSizeX(), seq.getSizeY()), canal);
        } else {
            return imageToArray(seq, Shape.make(seq.getSizeX(), seq.getSizeY(), seq.getSizeZ()), canal);
        }
    }

    public static ShapedArray imageToArray(Sequence seq, Shape shape, int canal) {
        if (seq.getSizeT() != 1) {
            throw new IllegalArgumentException("The input canno't be a 4D sequence");
        }
        try {
            double[] out = Array1DUtil.arrayToDoubleArray(seq.getDataCopyXYZT( canal ), seq.isSignedDataType());
            if (shape.rank() == 2) {    //2D input -> Rank = 2
                return Double2D.wrap(out, shape);
            } else {
                return Double3D.wrap(out, shape);
            }
        } catch (Exception e) {
            System.err.println("Canno't take only one canal "+e);
            if (shape.rank() == 2) {    //2D input -> Rank = 2
                return Double2D.wrap(icyImage3DToArray1D(seq.getAllImage(), shape.dimension(0), shape.dimension(1), 1, false), shape);
            } else {
                return Double3D.wrap(icyImage3DToArray1D(seq.getAllImage(), shape.dimension(0), shape.dimension(1), shape.dimension(2), false), shape);
            }
        }
    }


    public static ShapedArray sequenceToArray(IcyBufferedImage image, int width,int height) {
        //First we try if we can convert the data directly by using Icy Capacities
        int[] shape = new int[]{width, height};
        try {
            return Double2D.wrap(image.getDataCopyCXYAsDouble(), shape);
        } catch (Exception e) {
            return BufferedImageUtils.imageToArray(image);
        }
    }

    public static double[] shiftIcyPsf3DToArray1D(ArrayList<IcyBufferedImage>listPSF,int width, int height, int sizeZ,  boolean isComplex) {
        double[] out;
        if (isComplex) {
            out = new double[width*height*sizeZ*2];
        } else {
            out = new double[width*height*sizeZ];
        }
        double[] psfIn = icyImage3DToArray1D(listPSF, width, height, sizeZ, isComplex);
        if (psfIn.length != out.length) {
            System.err.println("Bad size for psf and output deconvutil l356");
        }
        CommonUtils.fftShift3D(psfIn,out, width, height, sizeZ);
        //CommonUtils.psf3DPadding1D(out, psfIn , width, height, sizeZ);
        return out;
    }

    public static ArrayList<IcyBufferedImage> arrayToImage(ShapedArray array) {
        Shape shape = array.getShape();
        int width = shape.dimension(0);
        int height = shape.dimension(1);
        if (array.getType() == Traits.DOUBLE) {
            if (shape.rank() == 2) {
                return fill(new IcyBufferedImage(width, height, ((Double2D)array).flatten()));
            } else if (shape.rank() == 3) {
                double[] data = ((Double3D)array).flatten();
                int sizeZ = shape.dimension(2);
                ArrayList<IcyBufferedImage> list = new ArrayList<IcyBufferedImage>();
                for (int j = 0; j < sizeZ; j++) {
                    double[] tmp = new double[width*height];
                    for (int i = 0; i < width*height; i++) {
                        tmp[i] = data[i+j*width*height];
                    }
                    IcyBufferedImage icyTmp = new IcyBufferedImage(width, height, tmp);
                    list.add(icyTmp);
                }
                return list;
            } else {
                throw new IllegalArgumentException("Rank of the Shaped Array can only be 2 or 3");
            }
        } else {
            if (shape.rank() == 2) {
                return fill(new IcyBufferedImage(width, height, ((Float2D)array).toFloat().flatten())); //Whatever the type other than float, we use float
            } else if (shape.rank() == 3) {
                float[] data = ((Float3D)array).flatten();
                int sizeZ = shape.dimension(2);
                ArrayList<IcyBufferedImage> list = new ArrayList<IcyBufferedImage>();
                for (int j = 0; j < sizeZ; j++) {
                    float[] tmp = new float[width*height];
                    for (int i = 0; i < width*height; i++) {
                        tmp[i] = data[i+j*width*height];
                    }
                    IcyBufferedImage icyTmp = new IcyBufferedImage(width, height, tmp);
                    list.add(icyTmp);
                }
                return list;
            } else {
                throw new IllegalArgumentException("Rank of the Shaped Array can only be 2 or 3");
            }
        }
    }

    private static ArrayList<IcyBufferedImage> fill(IcyBufferedImage img){
        ArrayList<IcyBufferedImage> list = new ArrayList<IcyBufferedImage>();
        list.add(img);
        return list;
    }

    /*
     * 
     * HERE BACKUP FROM DELETED FUNCTION FROM DECONVUTILS
     * 
     */

    //FIXME TMP COPY
    @Deprecated
    public static double[] icyImage3DToArray1D(ArrayList<IcyBufferedImage>listImage, int width, int height, int sizeZ, boolean isComplex) {
        double[] out;
        if (isComplex) {
            out = new double[2*sizeZ*width*height];
            int strideW = width;
            int strideH = width*height;
            for (int k = 0; k < sizeZ; k++) {
                double[] tmp = CommonUtils.imageToArray1D(listImage.get(k), false);
                for (int j = 0; j < height; j++) {
                    for (int i = 0; i < width; i++) {
                        out[2*i+2*j*strideW+2*k*strideH] = tmp[i+j*strideW];
                    }
                }
            }
        } else {
            out = new double[sizeZ*width*height];
            for (int j = 0; j < sizeZ; j++) {
                double[] tmp = CommonUtils.imageToArray1D(listImage.get(j), false);
                for (int i = 0; i < tmp.length; i++) {
                    out[i+j*tmp.length] = tmp[i];
                }
            }
        }
        return out;
    }

    public ArrayList<BufferedImage> arrayToIcyImage3D(double[] array, int job, boolean isComplex, int width, int height, int sizeZ){
        ArrayList<BufferedImage> out = new ArrayList<BufferedImage>();
        if (isComplex) {
            for (int k = 0; k < sizeZ; k++) {
                double[] tmp = new double[width*height];
                for (int j = 0; j < height; j++) {
                    for (int i = 0; i < width; i++) {
                        tmp[i+j*width] = array[2*i+2*j*width+2*k*height*width];
                    }
                }
                out.add(new IcyBufferedImage(width, height, tmp));
            }
        }else{
            for (int j = 0; j < sizeZ; j++) {
                double[] tmp = new double[width*height];
                for (int i = 0; i < width*height; i++) {
                    tmp[i] = array[i+j*height*width];
                }
                out.add(new IcyBufferedImage(width, height, tmp));
            }
        }
        //IMPORTANT WE DEPAD AS WE COMPUTE OR NOT ...
        //return CommonUtils.imageUnPad(out, sizePadding);
        return out;
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