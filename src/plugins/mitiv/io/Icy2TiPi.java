package plugins.mitiv.io;

import icy.image.IcyBufferedImage;
import icy.sequence.DimensionId;
import icy.sequence.Sequence;
import icy.type.DataType;
import mitiv.array.Array3D;
import mitiv.array.Array4D;
import mitiv.array.ArrayFactory;
import mitiv.array.ShapedArray;
import mitiv.base.Shape;

public class Icy2TiPi {
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
    public static ShapedArray sequenceToArray(Sequence seq,int c,int t) {
        return sequenceToArray( seq,c, -1, t);
    }
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
                data =  seq.getDataCopyXY(t, z, c); // FIXME it can be  seq.getDataXY(t, z, c);
                break;

            default:
                throw  new IllegalArgumentException("CZT Selection impossible");
        }

        return ArrayFactory.wrap(data, shape);
    }

    public static Sequence arrayToSequence(ShapedArray array)
    {
        return arrayToSequence( array,null);
    }

    public static Sequence arrayToSequence(ShapedArray array,Sequence sequence)
    {
        if (sequence == null )  {
            sequence = new Sequence();
        }

        switch (array.getRank()) {
            case 1:
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

}
