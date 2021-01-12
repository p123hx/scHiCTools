//
// Created by Bean Juice on 06/01/2021.
//

#include <xtensor/xsort.hpp>
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include "xtensor-blas/xlinalg.hpp"

using namespace std;

xt::xarray<double> PCA(xt::xarray<double> X, int dim = 2) {
    X = X - xt::mean(X, 0);
    xt::xarray<double> U, S, V;
    tie(U, S, V) = xt::linalg::svd(X, false, true);
    return xt::transpose(
            xt::linalg::dot(xt::view(V, xt::range(xt::placeholders::_, dim), xt::all()),
                            xt::transpose(X)));
}

//xt::xarray<double> MDS(xt::xarray<double> mat, int n = 2) {
//    xt::xarray<double> h = xt::eye(mat.size()) - xt::full_like(mat,1.);
//    xt::xarray <double> k = -0.5*xt::linalg::dot(xt::linalg::dot(h, mat*mat),h);
//    k[xt::isnan(k)] = 0;
//    xt::xarray <double> w,v;
//
//    tie(w,v) = xt::linalg::eig(k);
//    auto max_ = xt::view(xt::argsort(w),xt::range(xt::placeholders::_,-n-1,-1));
//    return xt::real(xt::col(v,max_));
//}

