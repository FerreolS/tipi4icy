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


import java.awt.image.BufferedImage;
import java.util.ArrayList;

import icy.image.IcyBufferedImage;
import icy.sequence.DimensionId;
import icy.sequence.Sequence;
import icy.type.collection.array.Array1DUtil;
import mitiv.array.Byte1D;
import mitiv.array.Byte2D;
import mitiv.array.Byte3D;
import mitiv.array.Byte4D;
import mitiv.array.Byte5D;
import mitiv.array.Double1D;
import mitiv.array.Double2D;
import mitiv.array.Double3D;
import mitiv.array.Double4D;
import mitiv.array.Double5D;
import mitiv.array.Float1D;
import mitiv.array.Float2D;
import mitiv.array.Float3D;
import mitiv.array.Float4D;
import mitiv.array.Float5D;
import mitiv.array.Int1D;
import mitiv.array.Int2D;
import mitiv.array.Int3D;
import mitiv.array.Int4D;
import mitiv.array.Int5D;
import mitiv.array.ShapedArray;
import mitiv.array.Short1D;
import mitiv.array.Short2D;
import mitiv.array.Short3D;
import mitiv.array.Short4D;
import mitiv.array.Short5D;
import mitiv.base.Shape;
import mitiv.base.Traits;
import mitiv.io.BufferedImageUtils;
import mitiv.utils.CommonUtils;

public class IcyBufferedImageUtils {
    static final byte czt = (byte)0b00000000;
    static final byte Czt = (byte)0b00000001;
    static final byte cZt = (byte)0b00000010;
    static final byte czT = (byte)0b00000100;
    static final byte cZT = (byte)0b00000110;
    static final byte CZt = (byte)0b00000011;
    static final byte CzT = (byte)0b00000101;
    static final byte CZT = (byte)0b00000111;
    /*Ferréol */

    public static ShapedArray sequenceToArray(Sequence seq) {
        return sequenceToArray( seq,-1, -1, -1);
    }
    public static ShapedArray sequenceToArray(Sequence seq,int c) {
        return sequenceToArray( seq,c, -1, -1);
    }
    public static ShapedArray sequenceToArray(Sequence seq,int c,int z) {
        return sequenceToArray( seq,c, z, -1);
    }
    public static ShapedArray sequenceToArray(Sequence seq,int c, int z, int t) {
        int nx, ny, nz, nc, nt;
        nx = seq.getSize(DimensionId.X);
        ny = seq.getSize(DimensionId.Y);
        nz = seq.getSize(DimensionId.Z);
        nc = seq.getSize(DimensionId.C);
        nt = seq.getSize(DimensionId.T);
        byte cztSelect = czt;
        int[] dims ={1,1,1,1,1};
        int ndims =0;
        if (c<0){
            if (nc>1){
                dims[ndims] = nc;
                ndims++;
            }
        }else{
            if (c<nc){
                cztSelect |= Czt;
            }else{
                throw  new IllegalArgumentException("Requested channel unavailable");
            }
        }
        if(nx>1){
            dims[ndims] = nx;
            ndims++;
        }
        if(ny>1){
            dims[ndims] = ny;
            ndims++;
        }
        if(nz>1){
            dims[ndims] = nz;
            ndims++;
        }
        if(nt>1){
            dims[ndims] = nt;
            ndims++;
        }
        // Remove singleton dimensions
        int[] newdims = new int[ndims];
        for (int i = 0; i < ndims; i++) {
            newdims[i] = dims[i];
        }

        Shape shape = new Shape(newdims);
        Object data;
        switch (cztSelect) {
            case czt:
                data =  seq.getDataCopyXYCZT();
                break;
            case Czt:
                data =  seq.getDataCopyXYZT(c);
                break;
            case cZt:
                throw  new IllegalArgumentException("seq.getDataCopyXYCT(z) is missing, ask Stephane for it");
                //                data =  seq.getDataCopyXYCT(z);
                //               break;
            case czT:
                data =  seq.getDataCopyXYCZ(t);
                break;
            case CZt:
                throw  new IllegalArgumentException("seq.getDataCopyXYT(c,z) is missing, ask Stephane for it");
                //                data =  seq.getDataCopyXYT(c,z);
                //               break;
            case CzT:
                data =  seq.getDataCopyXYZ(t, c);
                break;
            case cZT:
                data =  seq.getDataCopyXYC(t, c);
                break;
            case CZT:
                data =  seq.getDataCopyXY(t, z, c); // FIXME it can be  seq.getDataXY(t, z, c);
                break;

            default:
                throw  new IllegalArgumentException("CZT Selection impossible");
        };
        // Array1DUtil.arrayToArray(in, out, signed)
        switch (seq.getDataType_().getJavaType())
        {
            case BYTE:
                return wrapObject((byte[]) data, shape);
            case SHORT:
                return wrapObject((short[]) data, shape);
            case INT:
                return wrapObject((int[]) data, shape);
            case FLOAT:
                return wrapObject((float[]) data, shape);
            case DOUBLE:
                return wrapObject((double[]) data, shape);
            default:
                return null;
        }
    }



    private static ShapedArray wrapObject(double[] data, Shape shape) {
        switch (shape.rank()) {
            case 1:
                return Double1D.wrap(data, shape);
            case 2:
                return Double2D.wrap(data, shape);
            case 3:
                return Double3D.wrap(data, shape);
            case 4:
                return Double4D.wrap(data, shape);
            case 5:
                return Double5D.wrap(data, shape);
            default:
                throw new IllegalArgumentException("The input dimension must be of rank <6");
        }
    }

    private static ShapedArray wrapObject(float[] data, Shape shape) {
        switch (shape.rank()) {
            case 1:
                return Float1D.wrap(data, shape);
            case 2:
                return Float2D.wrap(data, shape);
            case 3:
                return Float3D.wrap(data, shape);
            case 4:
                return Float4D.wrap(data, shape);
            case 5:
                return Float5D.wrap(data, shape);
            default:
                throw new IllegalArgumentException("The input dimension must be of rank <6");
        }
    }

    private static ShapedArray wrapObject(int[] data, Shape shape) {
        switch (shape.rank()) {
            case 1:
                return Int1D.wrap(data, shape);
            case 2:
                return Int2D.wrap(data, shape);
            case 3:
                return Int3D.wrap(data, shape);
            case 4:
                return Int4D.wrap(data, shape);
            case 5:
                return Int5D.wrap(data, shape);
            default:
                throw new IllegalArgumentException("The input dimension must be of rank <6");
        }
    }


    private static ShapedArray wrapObject(short[] data, Shape shape) {
        switch (shape.rank()) {
            case 1:
                return Short1D.wrap(data, shape);
            case 2:
                return Short2D.wrap(data, shape);
            case 3:
                return Short3D.wrap(data, shape);
            case 4:
                return Short4D.wrap(data, shape);
            case 5:
                return Short5D.wrap(data, shape);
            default:
                throw new IllegalArgumentException("The input dimension must be of rank <6");
        }
    }

    private static ShapedArray wrapObject(byte[] data, Shape shape) {
        switch (shape.rank()) {
            case 1:
                return Byte1D.wrap(data, shape);
            case 2:
                return Byte2D.wrap(data, shape);
            case 3:
                return Byte3D.wrap(data, shape);
            case 4:
                return Byte4D.wrap(data, shape);
            case 5:
                return Byte5D.wrap(data, shape);
            default:
                throw new IllegalArgumentException("The input dimension must be of rank <6");
        }
    }


    /* Jonathan */
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
            return imageToArray(seq, new Shape(seq.getSizeX(), seq.getSizeY()), canal);
        } else {
            return imageToArray(seq, new Shape(seq.getSizeX(), seq.getSizeY(), seq.getSizeZ()), canal);
        }
    }

    public static ShapedArray imageToArray(Sequence seq, Shape shape, int canal) {
        if (seq.getSizeT() != 1) {
            throw new IllegalArgumentException("The input cannot be a 4D sequence");
        }
        try {
            double[] out = Array1DUtil.arrayToDoubleArray(seq.getDataCopyXYZT( canal ), seq.isSignedDataType());
            //Array3D. seq.getDataCopyXYZ(1, canal );
            /*   seq.getDataXYZ(1, canal);
            seq.getDataType_();
            Array2D ar;
            ar.wrap(seq.getDataCopyXYZT( canal ), shape);*/
            // ShapedArray
            if (shape.rank() == 2) {    //2D input -> Rank = 2
                return Double2D.wrap(out, shape);
            } else {
                return Double3D.wrap(out, shape);
            }
        } catch (Exception e) {
            System.err.println("Cannot take only one canal "+e);
            if (shape.rank() == 2) {    //2D input -> Rank = 2
                return Double2D.wrap(icyImage3DToArray1D(seq.getAllImage(), shape.dimension(0), shape.dimension(1), 1, false), shape);
            } else {
                return Double3D.wrap(icyImage3DToArray1D(seq.getAllImage(), shape.dimension(0), shape.dimension(1), shape.dimension(2), false), shape);
            }
        }
    }


    /*  public static ShapedArray sequenceToArray(IcyBufferedImage image, int width,int height) {
        //First we try if we can convert the data directly by using Icy Capacities
        int[] shape = new int[]{width, height};
        try {
            return Double2D.wrap(image.getDataCopyCXYAsDouble(), shape);
        } catch (Exception e) {
            return BufferedImageUtils.imageToArray(image);
        }
    }
     */
    public static Sequence IcyBufferedToSequence(ArrayList<IcyBufferedImage> image) {
        Sequence tmp = new Sequence();
        for (int i = 0; i < image.size(); i++) {
            tmp.addImage(i, image.get(i));
        }
        return tmp;
    }

    public static Sequence BufferedToSequence(ArrayList<BufferedImage> image) {
        Sequence tmp = new Sequence();
        for (int i = 0; i < image.size(); i++) {
            tmp.addImage(i, image.get(i));
        }
        return tmp;
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

    public static ArrayList<BufferedImage> arrayToIcyImage3D(double[] array, boolean isComplex, int width, int height, int sizeZ){
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

    public static Sequence arrayToSequence(double[] array, boolean isComplex, int width, int height, int sizeZ){
        ArrayList<BufferedImage> tmp = arrayToIcyImage3D(array, isComplex, width, height, sizeZ);
        return BufferedToSequence(tmp);
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