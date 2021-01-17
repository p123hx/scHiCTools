//
// Created by Bean Juice on 08/01/2021.
//

#include "processing_utils.h"

using namespace std;


xt::xarray<double> matrix_operation(xt::xarray<double> mat, vector<string> operations) {
    xt::xarray<double> tmp;
    for (string op : operations) {
        transform(op.begin(), op.end(), op.begin(), ::tolower);
        if (op == "convolution") {
            tmp = convolution(mat, 3);
        } else {
            throw "Not implemented";
        }

    }
    return tmp;
}

xt::xarray<double> convolution(xt::xarray<double> mat, int kernel_shape) {
//credit: https://github.com/chaowang15/fast-image-convolution-cpp/blob/master/src/convolution.cpp

    xt::xarray<double> conv = xt::ones<double>({kernel_shape, kernel_shape}) / pow
            (kernel_shape, 2.0);

    xt::xarray<double> outMat = xt::zeros_like(mat);

    double temp;
    int k, l;
    int low_k, high_k, low_l, high_l;

    int h_dst = mat.shape(0), w_dst = mat.shape(1);
    int h_src = h_dst, w_src = w_dst;
    for (int i = 0; i < h_dst; ++i) {
        low_k = max(0, i - int((kernel_shape - 1.0) / 2.0));
        high_k = min(h_src - i, i + int(kernel_shape / 2.0));
        for (int j = 0; j < w_dst; j++) {
            low_l = std::max(0, j - int((kernel_shape - 1.0) / 2.0));
            high_l = std::min(w_src - 1, j + int(kernel_shape / 2.0));
            temp = 0.0;
            for (k = low_k; k <= high_k; ++k) {
                for (l = low_l; l <= high_l; ++l) {
                    temp += mat(k * w_src + l) *
                            conv((i - k + int(kernel_shape / 2.0)) *
                                 kernel_shape +
                                 (j - l + int(kernel_shape / 2.0)));
                }
            }
            outMat(i * w_dst + j) = temp;
        }
    }


//    int i, j, m, n;
//    int inInd, inInd2, outInd = 0, kInd = 0;
//    int kCenterX = kernel_shape >> 1, kCenterY = kernel_shape >> 1;
//    int rowMin, rowMax;  // to check boundary of input array
//    int colMin, colMax;  //
//    int dataSizeX, dataSizeY = mat.dimension(); wrong!
//    inInd = inInd2 = dataSizeX * kCenterY + kCenterX;
//    for (i = 0; i < dataSizeY; ++i) {
//        rowMax = i + kCenterY;
//        rowMin = i - dataSizeY + kCenterY;
//        for (j = 0; j < dataSizeX; ++j) {
//            colMax = j + kCenterX;
//            colMin = j - dataSizeX + kCenterX;
//
//            for (m = 0; m < kernel_shape; ++m) {
//                if (m <= rowMax && m > rowMin) {
//                    for (n = 0; n < kernel_shape; ++n) {
//                        // check the boundary of array
//                        if (n <= colMax && n > colMin)
//                            outMat(outInd) += mat(inInd - n) * conv(kInd);
//                        ++kInd;  // next kernel
//                    }
//                } else {
//                    kInd += kernel_shape;
//                }
//                inInd -= dataSizeX;
//            }
//            kInd = 0;
//            inInd = ++inInd2;
//            ++outInd;
//        }
//    }
    return outMat;

}