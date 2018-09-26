/*  Copyright (C) 2017  Ferreol Soulez ferreol.soulez@univ-lyon1.fr
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  */

package plugins.mitiv.io;

import icy.image.IcyBufferedImage;
import icy.sequence.DimensionId;
import icy.sequence.Sequence;
import icy.type.DataType;
import mitiv.array.Array2D;
import mitiv.array.Array3D;
import mitiv.array.Array4D;
import mitiv.array.Array5D;
import mitiv.array.ArrayFactory;
import mitiv.array.ShapedArray;
import mitiv.base.Shape;

/**
 * TiPi ShapedArray <->Icy Sequence conversion class
 * @author Ferréol ferreol.soulez@epfl.ch
 * @see mitiv.array.ShapedArray
 * @see icy.sequence
 *
 */
public class Icy2TiPi {
    static final byte czt = (byte) 0x0; // 0b00000000;
    static final byte Czt = (byte) 0x1; // 0b00000001;
    static final byte cZt = (byte) 0x2; // 0b00000010;
    static final byte czT = (byte) 0x4; // 0b00000100;
    static final byte cZT = (byte) 0x6; // 0b00000110;
    static final byte CZt = (byte) 0x3; // 0b00000011;
    static final byte CzT = (byte) 0x5; // 0b00000101;
    static final byte CZT = (byte) 0x7; // 0b00000111;
    /*Ferréol */


    /**
     * Convert sequence to array
     * @param seq input sequence
     * @return ShapedArray
     */
    public static ShapedArray sequenceToArray(Sequence seq) {
        return sequenceToArray( seq,-1, -1, -1);
    }
    /**
     * Convert sequence to array
     * @param   seq     input sequence
     * @param   c       channel index
     * @return  ShapedArray
     */
    public static ShapedArray sequenceToArray(Sequence seq,int c) {
        return sequenceToArray( seq,c, -1, -1);
    }
    /**
     * Convert sequence to array
     * @param seq   input sequence
     * @param c     channel index
     * @param t     time index
     * @return ShapedArray
     */
    public static ShapedArray sequenceToArray(Sequence seq,int c,int t) {
        return sequenceToArray( seq,c, -1, t);
    }
    /**
     * Convert sequence to array
     * @param seq input sequence
     * @param c     channel index
     * @param z     depth index
     * @param t     time index
     * @return ShapedArray
     */
    public static ShapedArray sequenceToArray(Sequence seq,int c, int z, int t) {
        int nx, ny, nz, nc, nt;


        switch (seq.getDataType_())
        {
            case BYTE:
                seq = icy.sequence.SequenceUtil.convertToType(seq, DataType.SHORT, false);
                break;
            case USHORT:
                seq = icy.sequence.SequenceUtil.convertToType(seq, DataType.INT, false);
                break;
            case UINT:
                seq = icy.sequence.SequenceUtil.convertToType(seq, DataType.LONG, false);
                break;
            case ULONG:
                seq = icy.sequence.SequenceUtil.convertToType(seq, DataType.DOUBLE, false);
                break;
            default:
                break;
        }

        /* extract sequence dimension */
        nx = seq.getSize(DimensionId.X);
        ny = seq.getSize(DimensionId.Y);
        nz = seq.getSize(DimensionId.Z);
        nc = seq.getSize(DimensionId.C);
        nt = seq.getSize(DimensionId.T);

        byte cztSelect = czt; // CZT selector

        int[] dims ={1,1,1,1,1};
        int ndims =0;
        if (c<0){ // No C selection
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

        if(z<0){ // no Z selection
            if(nz>1){
                dims[ndims] = nz;
                ndims++;
            }
        }else{
            if (z<nz){
                cztSelect |= cZt;
            }else{
                throw  new IllegalArgumentException("Requested z unavailable");
            }
        }

        if(t<0){ // no T selection
            if(nt>1){
                dims[ndims] = nt;
                ndims++;
            }
        }else{
            if (t<nt){
                cztSelect |= czT;
            }else{
                throw  new IllegalArgumentException("Requested z unavailable");
            }
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
                data =  seq.getDataCopyXYC(t, z);
                break;
            case CZT:
                data =  seq.getDataCopyXY(t, z, c); //  it can be  seq.getDataXY(t, z, c);
                break;

            default:
                throw  new IllegalArgumentException("CZT Selection impossible");
        }

        return ArrayFactory.wrap(data, shape);
    }

    /**
     * Convert ShapedArray to Sequence
     * @param array     ShapedArray
     * @return          Sequence
     */
    public static Sequence arrayToSequence(ShapedArray array)
    {
        return arrayToSequence( array,null);
    }

    /**
     * Copy a ShapedArray into a Sequence
     * Create a new sequence if the input sequence is null
     * @param array     ShapedArray
     * @param sequence  Input sequence
     * @return          Resulting sequence
     */
    public static Sequence arrayToSequence(ShapedArray array,Sequence sequence)
    {
        if (sequence == null )  {
            sequence = new Sequence();
        }
        switch (array.getRank()) {
            case 1:
                sequence.setImage(0,0, new IcyBufferedImage(array.getDimension(0), 1, array.flatten(),true,false));

                sequence.setImage(0,0, new IcyBufferedImage(array.getDimension(0), 1, array.flatten(),true,false));
                break;
            case 2:
                sequence.setImage(0,0, new IcyBufferedImage(array.getDimension(0), array.getDimension(1), array.flatten(),true,false));
                break;
            case 3:
                for (int j = 0; j < array.getDimension(2); j++) {
                    sequence.setImage(0,j, new IcyBufferedImage(array.getDimension(0), array.getDimension(1),((Array3D)array).slice(j).flatten() ,true,false));
                }
                break;

            case 4:

                for (int k = 0; k < array.getDimension(3); k++) {
                    for (int j = 0; j < array.getDimension(2); j++) {
                        sequence.setImage(k,j, new IcyBufferedImage(array.getDimension(0), array.getDimension(1),((Array4D)array).slice(k).slice(j).flatten() ,true,false));
                    }
                }
            default:
                throw new IllegalArgumentException(" arrayToSequence can convert only 1D to 4D arrays");
        }
        return sequence;
    }


    /**
     * Copy a ShapedArray into a Sequence
     * Create a new sequence if the input sequence is null
     * @param array     ShapedArray
     * @param channelIndex  dimension that corresponds to the channel
     * @param sequence  Input sequence
     * @return          Resulting sequence
     */
    public static Sequence arrayToSequence(ShapedArray array,int channelIndex, Sequence sequence)
    {

        if (channelIndex <0){
            return arrayToSequence( array,sequence);
        }else if (channelIndex>array.getRank()){
            throw new IllegalArgumentException(" The channel index cannot be larger than the array rank");
        }

        if (sequence == null )  {
            sequence = new Sequence();
        }

        int newdims[] = new int[array.getRank()-1];
        {
            int k=0,n=0;
            while (k<array.getRank()-1) {
                if (n!=channelIndex){
                    newdims[k] =    array.getDimension(n);
                    k++;
                }
                n++;
            }
        }

        double[][] t = new double[array.getDimension(channelIndex)][newdims[0]* newdims[1]];

        switch (array.getRank()-1) {
            case 1:
                for(int n=0;n<array.getDimension(channelIndex);++n){
                    t[n] =  ((Array2D) array).slice(n,channelIndex).toDouble().flatten();
                }
                sequence.setImage(0,0, new IcyBufferedImage(newdims[0], 1, t,true,false));
                break;
            case 2:
                for(int n=0;n<array.getDimension(channelIndex);++n){
                    t[n] =  ((Array3D) array).slice(n,channelIndex).toDouble().flatten();
                }

                sequence.setImage(0,0, new IcyBufferedImage(newdims[0], newdims[1],t,true,false));
                break;
            case 3:
                for (int j = 0; j < newdims[2]; j++) {
                    for(int n=0;n<array.getDimension(channelIndex);++n){
                        t[n] =  ((Array4D) array).slice(n,channelIndex).slice(j).toDouble().flatten();
                    }
                    sequence.setImage(0,j, new IcyBufferedImage(newdims[0], newdims[1],t ,true,false));
                }
                break;

            case 4:

                for (int k = 0; k < array.getDimension(3); k++) {
                    for (int j = 0; j < array.getDimension(2); j++) {for(int n=0;n<array.getDimension(channelIndex);++n){
                        t[n] =  ((Array5D) array).slice(n,channelIndex).slice(k).slice(j).toDouble().flatten();
                    }
                    sequence.setImage(k,j, new IcyBufferedImage(newdims[0], newdims[1],t ,true,false));
                    }
                }
            default:
                throw new IllegalArgumentException(" arrayToSequence can convert only 1D to 4D arrays");
        }
        return sequence;
    }


}
