
//
// Created by Bean Juice on 16/12/2020.
//

#include "reproducibility.h"

#include <utility>
#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xadapt.hpp"
#include "xtensor-blas/xlinalg.hpp"


using namespace std;
using namespace std::chrono;

double pearsoncoeff(xt::xarray<double> x, xt::xarray<double> y, double size) {
    return xt::sum((x - xt::mean(x)) * (y - xt::mean(y)))() /
           (size * xt::stddev(x)() * xt::stddev(y)());
}

xt::xarray<double> euc_pdist_square(xt::xarray<double> x) {
    int rol, col = x.dimension();
    vector<double> tmp(rol * (rol - 1) / 2);
    for (int i = 0; i < rol; i++) {
        for (int j = 0; j < i; j++) {
            double dum = 0;
            for (int k = 0; k < col; k++) {
                dum += sqrt(pow(x(i, k) - x(j, k), 2));
            }
            tmp.push_back(dum);
        }
    }
    xt::xarray<double> ans = xt::zeros<double>({rol, rol});

    for (int i = 0; i < rol; i++)
        for (int j = 0; j < i; j++) {
            int ind = rol * (rol - 1) / 2 - (rol - i) * (rol - i - 1) / 2 + j - i - 1;
            ans(j, i) = tmp[ind];
            ans(i, j) = tmp[ind];
        }
    return ans;
}

//Note that this selfdefined zscore func is using propagating policy on axis
xt::xarray<double> zscore_prop(xt::xarray<double> a, int axis) {
    xt::xarray<double> mns = mean(a, {axis});
    xt::xarray<double> sstd = stddev(a, {axis});
    return (a - mns) / sstd;
}

xt::xarray<double>
pairwise_distance(xt::xarray<xt::xarray<double>> all_strata, string similarity_method,
                  bool print_time = false, double sigma = .5,
                  unsigned window_size = 10) {
    transform(similarity_method.begin(), similarity_method.end(),
              similarity_method.begin(), ::tolower);
    auto t0 = high_resolution_clock::now();
    high_resolution_clock::time_point t1, t2;
    xt::xarray<double> zscores;

    int n_cells, n_bins = all_strata[0].dimension();
    int all_size = all_strata.size();
    xt::xarray<xt::xarray<double>> tmp = xt::zeros<double>({1, all_size});
    xt::xarray<double> distance_mat;

    if (similarity_method == "inner_product" or similarity_method == "innerproduct") {
        int index = 0;
        for (xt::xarray<double> stratum : all_strata) {
            xt::xarray<double> z = zscore_prop(stratum, 1);
            //??
            z[isnan(z)] = 0;
            tmp(index) = z;
            index++;
        }
        zscores = xt::concatenate(xt::xtuple(tmp), 1);
        t1 = high_resolution_clock::now();
        xt::xarray<double> inner = xt::linalg::dot(zscores, xt::transpose(zscores));
        inner[inner > 1] = 1;
        inner[inner < -1] = -1;
        xt::xarray<double> distance_mat = xt::sqrt(2 - 2 * inner);
        t2 = high_resolution_clock::now();
    } else if (similarity_method == "hicrep") {

        int n_strata = all_strata.size();
        xt::xarray<double> weighted_std = xt::zeros<double>({n_cells, n_bins});
        int i = 0;
        for (auto stratum : all_strata) {
            xt::xarray<double> mean = xt::mean(stratum, 1)();
            xt::xarray<double> std = xt::stddev(stratum)();
            xt::col(weighted_std, i) = sqrt(n_bins - i) * std;
            all_strata[i] -= xt::col(mean, NULL);
            i++;
        }
        xt::xarray<double> scores = xt::concatenate(xt::xtuple(all_strata), 1);
        t1 = high_resolution_clock::now();
        xt::xarray<double> inner = xt::linalg::dot(scores, xt::transpose(scores));
        inner[inner > 1] = 1;
        inner[inner < -1] = -1;
        distance_mat = xt::sqrt(2 - 2 * inner);
        t2 = high_resolution_clock::now();

    } else if (similarity_method == "old_hicrep") {

        xt::xarray<double> similarity = xt::ones<double>({n_cells, n_bins});
        int ind_old=0;
        for (int i = 0; i < n_cells; i++) {
            for (int j = i + 1; j < n_cells; j++) {
                xt::xarray<double> corrs, weights = xt::empty<double>({1,all_size});
                for (auto stratum : all_strata) {

                    xt::xarray<double> s1 = xt::row(stratum, i);
                    xt::xarray<double> s2 = xt::row(stratum, j);
                    if (xt::variance(s1)() == 0 or xt::variance(s2)() == 0) {
                        weights(ind_old)=0;
                        corrs(ind_old)=0;
                    } else {
                        weights(ind_old)=(n_bins * xt::stddev(s1)() * xt::stddev(s2)());
                        corrs(ind_old)=(pearsoncoeff(s1, s2, n_bins));

                    }
                    ind_old++;
                }
                corrs = xt::nan_to_num(corrs);
                double s = xt::linalg::dot(corrs, weights)() / (xt::sum(weights)());
                similarity[i, j] = s;
                similarity[j, i] = s;
            }
        }
        t1 = high_resolution_clock::now();
        distance_mat = xt::sqrt(2 - 2 * similarity);
        t2 = high_resolution_clock::now();
    } else if (similarity_method == "selfish") {
        int n_windows = n_bins / window_size;
        xt::xarray<double> all_windows = xt::zeros<double>({n_cells, n_bins});
        int i = 0;
        for (auto stratum : all_strata) {
            for (int j = 0; j < n_windows; j++) {

                xt::col(all_windows, j) += xt::sum(xt::view(stratum, xt::all(),
                                                            xt::range(j * window_size,
                                                                      (j + 1) *
                                                                      window_size - i)),
                                                   1);
            }
        }
        t1 = high_resolution_clock::now();

        xt::xarray<double> fingerprints = xt::zeros<double>(
                {n_cells, n_windows * (n_windows - 1) / 2});
        int k = 0;
        for (int i = 0; i < n_windows; i++)
            for (int j = 0; j < n_windows - i - 1; j++) {
                xt::col(fingerprints, k) =
                        xt::col(all_windows, i) > xt::col(all_windows, j);
                k += 1;
            }
        xt::xarray<double> distance = euc_pdist_square(fingerprints);

        distance_mat = xt::sqrt(2 - 2 * xt::exp(-sigma * distance));
        t2 = high_resolution_clock::now();
    } else {
        throw "Method {0} not supported. Only \"inner_product\", \"HiCRep\", \"old_hicrep\" and \"Selfish\".";
    }


    auto duration1 = duration_cast<microseconds>(t1 - t0);
    auto duration2 = duration_cast<microseconds>(t2 - t1);
    if (print_time) {
        cout << "Time 1:" << duration1.count() << endl
             << "Time 2:" << duration2.count() << endl;
    }

    return distance_mat;
}


