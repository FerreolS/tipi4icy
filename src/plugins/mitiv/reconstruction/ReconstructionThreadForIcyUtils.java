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

import icy.image.IcyBufferedImage;
import icy.type.DataType;
import mitiv.array.ByteArray;
import mitiv.array.DoubleArray;
import mitiv.array.FloatArray;
import mitiv.array.IntArray;
import mitiv.array.ShapedArray;
import mitiv.array.ShortArray;
import mitiv.base.Traits;

public class ReconstructionThreadForIcyUtils {

    public static DataType toIcyDataType(int type) {
        switch (type) {
        case Traits.BYTE: return DataType.BYTE;
        case Traits.SHORT: return DataType.SHORT;
        case Traits.INT: return DataType.INT;
        case Traits.LONG: return DataType.LONG;
        case Traits.FLOAT: return DataType.FLOAT;
        case Traits.DOUBLE: return DataType.DOUBLE;
        }
        throw new IllegalArgumentException("Unknown data type.");
    }

    public void fetchResult(ShapedArray result, IcyBufferedImage image) {
        int n = result.getNumber();
        switch (result.getType()) {
        case Traits.BYTE:
        {
            byte[] src = ((ByteArray)result).flatten(false);
            byte[] dst = image.getDataXYAsByte(0);
            for (int j = 0; j < n; ++j) {
                dst[j] = src[j];
            }
            break;
        }
        case Traits.SHORT:
        {
            short[] src = ((ShortArray)result).flatten(false);
            short[] dst = image.getDataXYAsShort(0);
            for (int j = 0; j < n; ++j) {
                dst[j] = src[j];
            }
            break;
        }
        case Traits.INT:
        {
            int[] src = ((IntArray)result).flatten(false);
            int[] dst = image.getDataXYAsInt(0);
            for (int j = 0; j < n; ++j) {
                dst[j] = src[j];
            }
            break;
        }
        case Traits.LONG:
            throw new IllegalArgumentException("Type 'long' not supported in Icy.");
        case Traits.FLOAT:
        {
            float[] src = ((FloatArray)result).flatten(false);
            float[] dst = image.getDataXYAsFloat(0);
            for (int j = 0; j < n; ++j) {
                dst[j] = src[j];
            }
            break;
        }
        case Traits.DOUBLE:
        {
            double[] src = ((DoubleArray)result).flatten(false);
            double[] dst = image.getDataXYAsDouble(0);
            for (int j = 0; j < n; ++j) {
                dst[j] = src[j];
            }
            break;
        }
        default:
            throw new IllegalArgumentException("Unknown type.");
        }
        image.dataChanged();
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
