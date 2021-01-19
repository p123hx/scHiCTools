//
// Created by Bean Juice on 16/12/2020.
//

#ifndef SCHICTOOLS_REPRODUCIBILITY_H
#define SCHICTOOLS_REPRODUCIBILITY_H

#include <utility>
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include "xtensor-blas/xlinalg.hpp"

using namespace std;
using namespace std::chrono;

double pearsoncoeff(xt::xarray<double> x, xt::xarray<double> y, double size);

xt::xarray<double> euc_pdist_square(xt::xarray<double> x);

xt::xarray<double> zscore_prop(xt::xarray<double> a, int axis);

pair<xt::xarray<double>,double>
pairwise_distance(vector<xt::xarray<double>> all_strata, string similarity_method,
                  bool print_time = false, double sigma = .5,
                  unsigned window_size = 10);

#endif //SCHICTOOLS_REPRODUCIBILITY_H
