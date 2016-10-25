package plugins.mitiv.io;

import icy.image.IcyBufferedImage;
import icy.sequence.DimensionId;
import icy.sequence.Sequence;
import icy.type.DataType;
import mitiv.array.Array3D;
import mitiv.array.Array4D;
import mitiv.array.ArrayFactory;
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
import mitiv.array.Long1D;
import mitiv.array.Long2D;
import mitiv.array.Long3D;
import mitiv.array.Long4D;
import mitiv.array.Long5D;
import mitiv.array.ShapedArray;
import mitiv.array.Short1D;
import mitiv.array.Short2D;
import mitiv.array.Short3D;
import mitiv.array.Short4D;
import mitiv.array.Short5D;
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
        };

        switch (seq.getDataType_().getJavaType())
        {
            case BYTE:
                return ArrayFactory.wrap((byte[]) data, shape);
            case SHORT:
                return ArrayFactory.wrap((short[]) data, shape);
            case INT:
                return ArrayFactory.wrap((int[]) data, shape);
            case LONG:
                return ArrayFactory.wrap((long[]) data, shape);
            case FLOAT:
                return ArrayFactory.wrap((float[]) data, shape);
            case DOUBLE:
                return ArrayFactory.wrap((double[]) data, shape);
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

    private static ShapedArray wrapObject(long[] data, Shape shape) {
        switch (shape.rank()) {
            case 1:
                return Long1D.wrap(data, shape);
            case 2:
                return Long2D.wrap(data, shape);
            case 3:
                return Long3D.wrap(data, shape);
            case 4:
                return Long4D.wrap(data, shape);
            case 5:
                return Long5D.wrap(data, shape);
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
