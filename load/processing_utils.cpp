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
    cout<<tmp<<endl;
    return tmp;
}

xt::xarray<double> convolution(xt::xarray<double> mat, int kernel_shape) {
    //see https://github.com/chaowang15/fast-image-convolution-cpp.git
    //TODO: check up why convolution ends up with all zeroes and make it right
    xt::xarray<double> conv = xt::ones<double>({kernel_shape, kernel_shape}) / pow
            (kernel_shape, 2.0);

    xt::xarray<double> outMat = xt::zeros_like(mat);
    int i, j, m, n;
    int inInd, inInd2, outInd = 0, kInd = 0;
    int kCenterX = kernel_shape >> 1, kCenterY = kernel_shape >> 1;
    int rowMin, rowMax;  // to check boundary of input array
    int colMin, colMax;  //
    int dataSizeX, dataSizeY = mat.dimension();
    inInd = inInd2 = dataSizeX * kCenterY + kCenterX; 
    for (i = 0; i < dataSizeY; ++i) {
        rowMax = i + kCenterY;
        rowMin = i - dataSizeY + kCenterY;
        for (j = 0; j < dataSizeX; ++j) {
            colMax = j + kCenterX;
            colMin = j - dataSizeX + kCenterX;

            for (m = 0; m < kernel_shape; ++m) {
                if (m <= rowMax && m > rowMin) {
                    for (n = 0; n < kernel_shape; ++n) {
                        // check the boundary of array
                        if (n <= colMax && n > colMin)
                            outMat(outInd) += mat(inInd - n) * conv(kInd);
                        ++kInd;  // next kernel
                    }
                } else {
                    kInd += kernel_shape;
                }
                inInd -= dataSizeX;
            }
            kInd = 0;
            inInd = ++inInd2;
            ++outInd;
        }
    }
    return outMat;

}