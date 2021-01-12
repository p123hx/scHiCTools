//
// Created by Bean Juice on 08/01/2021.
//

#ifndef SCHICTOOLS_PROCESSING_UTILS_H
#define SCHICTOOLS_PROCESSING_UTILS_H
#include <cmath>
#include <set>
#include <map>
#include <string>
#include <fstream>
#include <iostream>
#include <list>
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include "xtensor-blas/xlinalg.hpp"

using namespace std;
xt::xarray<double> matrix_operation(xt::xarray<double>mat,vector<string>operations);
xt::xarray<double>convolution(xt::xarray<double>mat,int kernel_shape =3);
#endif //SCHICTOOLS_PROCESSING_UTILS_H
