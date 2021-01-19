
//
// Created by Hongxi on 16/12/2020.
//
#include <utility>
#include <fstream>
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include "xtensor-blas/xlinalg.hpp"
#include "float.h"

using namespace std;
using namespace std::chrono;

double pearsoncoeff(xt::xarray<double> x, xt::xarray<double> y, double size) {
    return xt::sum((x - xt::mean(x)) * (y - xt::mean(y)))() /
           (size * xt::stddev(x)() * xt::stddev(y)());
}

//along axis = 1
xt::xarray<double> concatenate_axis1(vector<xt::xarray<double>> all_strata) {
    if (all_strata.empty()) { throw "all_strata is an empty vector. Contatenation fail"; }
    xt::xarray<double> ans = all_strata[0];
    for (int i = 1; i < all_strata.size(); i++) {
        ans = xt::concatenate(xt::xtuple(ans, all_strata[i]), 1);
    }
    return ans;
}

xt::xarray<double> euc_pdist_square(xt::xarray<double> x) {
    int row = x.shape(0), col = x.shape(1);
    vector<double> tmp(row * (row - 1) / 2);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < i; j++) {
            double dum = 0;
            for (int k = 0; k < col; k++) {
                dum += sqrt(pow(x(i, k) - x(j, k), 2));
            }
            tmp.push_back(dum);
        }
    }
    xt::xarray<double> ans = xt::zeros<double>({row, row});

    for (int i = 0; i < row; i++)
        for (int j = 0; j < i; j++) {
            int ind = row * (row - 1) / 2 - (row - i) * (row - i - 1) / 2 + j - i - 1;
            ans(j, i) = tmp[ind];
            ans(i, j) = tmp[ind];
        }
    return ans;
}

//Note that this selfdefined zscore func is using propagating policy on axis
xt::xarray<double> zscore_prop(xt::xarray<double> a, int axis) {
    xt::xarray<double> mns = xt::mean(a, {axis});
    cout<<"Mean:\n"<<mns;
    xt::xarray<double> sstd = xt::stddev(a, {axis});
    cout<<"\nstd:\n"<<sstd;
    int row = a.shape(0),col=a.shape(1);
    for(int i=0;i<row;i++){
        for(int j=0;j<col;j++){
            a(i,j)=(a(i,j)-mns(i))/sstd(i);
        }
    }
    return a;
}

pair<xt::xarray<double>,double>
pairwise_distance(vector<xt::xarray<double>> all_strata, string similarity_method,
                  bool print_time = true, double sigma = .5,
                  unsigned window_size = 10) {
    transform(similarity_method.begin(), similarity_method.end(),
              similarity_method.begin(), ::tolower);
    auto t0 = high_resolution_clock::now();
    high_resolution_clock::time_point t1, t2;
    xt::xarray<double> zscores;

    int n_cells = all_strata[0].shape(0), n_bins = all_strata[0].shape(1);
    int all_size = all_strata.size();
    vector<xt::xarray<double>> tmp;
    xt::xarray<double> distance_mat;

    if (similarity_method == "inner_product" or similarity_method == "innerproduct") {

        for (xt::xarray<double> stratum : all_strata) {
            cout<<"stratum: \n"<<stratum;
            cout<<stratum.shape(0)<<" by "<<stratum.shape(1)<<endl;
            xt::xarray<double> z = zscore_prop(stratum, 1);
            //??
            cout<<"\nz:\n"<<z;
            for (int i = 0; i < z.size(); ++i) {
                if(isnan(z(i))) z(i)=0.0;
            }
            tmp.push_back(z);

        }
        zscores = concatenate_axis1(tmp);
        cout<<"zscores:\n"<<zscores;
        t1 = high_resolution_clock::now();
        xt::xarray<double> inner = xt::linalg::dot(zscores, xt::transpose(zscores))
                /zscores.shape(1);
        for (int i = 0; i < inner.size(); i++) {
            if (inner(i) > 1) inner(i) = 1.0;
            if (inner(i) < -1) inner(i) = -1.0;
        }
        distance_mat = xt::sqrt(2 - 2 * inner);
        t2 = high_resolution_clock::now();
    } else if (similarity_method == "hicrep") {

        int n_strata = all_strata.size();
        xt::xarray<double> weighted_std = xt::zeros<double>({n_cells, n_bins});
        int i = 0;
        for (auto stratum : all_strata) {
//            cout << "i=" << i << endl;
            xt::xarray<double> mean = xt::mean(stratum, {1});
            xt::xarray<double> std = xt::stddev(stratum, {1});
//            cout << "Mean:\n" << mean << endl << "std:\n" << std << "\nstrata "
//                                                                    "be4\n" << stratum
//                 << "\nstrata after\n";
            xt::col(weighted_std, i) = sqrt(n_bins - i) * std;
            for (int j = 0; j < mean.size(); j++) {
                xt::row(all_strata[i], j) = xt::row(stratum, j) - mean(j);
            }

//            cout << all_strata[i] << endl;
            i++;
        }
        xt::xarray<double> scores = concatenate_axis1(all_strata);
//        cout << "scores:\n " << scores << endl;
//        cout << "weighted_std:\n" << weighted_std << endl;
        t1 = high_resolution_clock::now();
        xt::xarray<double> inner = xt::linalg::dot(scores, xt::transpose(scores)) /
                                   (xt::linalg::dot(weighted_std,
                                                    xt::transpose(weighted_std)) +
                                    DBL_MIN);
//        cout << "dinominator:\n" << xt::linalg::dot(scores, xt::transpose(scores))
//             << "\nnumeritor\n"
//             << (xt::linalg::dot(weighted_std, xt::transpose(weighted_std)) + DBL_MIN
//             ) << "\n inner be4\n" << inner;
        for (int i = 0; i < inner.size(); i++) {
            if (inner(i) > 1) inner(i) = 1.0;
            if (inner(i) < -1) inner(i) = -1.0;
        }
//        cout << "\ninner: " << inner << endl;
        distance_mat = xt::sqrt(2 - 2 * inner);
//        cout << "distance_mat infunc\n" << distance_mat << endl;
        t2 = high_resolution_clock::now();

    } else if (similarity_method == "old_hicrep") {

        xt::xarray<double> similarity = xt::ones<double>({n_cells, n_bins});
        int ind_old = 0;
        for (int i = 0; i < n_cells; i++) {
            for (int j = i + 1; j < n_cells; j++) {
                xt::xarray<double> corrs, weights = xt::empty<double>({1, all_size});
                for (auto stratum : all_strata) {

                    xt::xarray<double> s1 = xt::row(stratum, i);
                    xt::xarray<double> s2 = xt::row(stratum, j);
                    if (xt::variance(s1)() == 0 or xt::variance(s2)() == 0) {
                        weights(ind_old) = 0;
                        corrs(ind_old) = 0;
                    } else {
                        weights(ind_old) = (n_bins * xt::stddev(s1)() *
                                            xt::stddev(s2)());
                        corrs(ind_old) = (pearsoncoeff(s1, s2, n_bins));

                    }
                    ind_old++;
                }
                corrs = xt::nan_to_num(corrs);
                double s = xt::linalg::dot(corrs, weights)() / (xt::sum(weights)());
                similarity(i, j) = s;
                similarity(i, j) = s;
            }
        }
        t1 = high_resolution_clock::now();
        distance_mat = xt::sqrt(2 - 2 * similarity);
        t2 = high_resolution_clock::now();
    } else if (similarity_method == "selfish") {
        int n_windows = n_bins / window_size;
        xt::xarray<double> all_windows = xt::zeros<double>({n_cells, n_bins});
        int stra_i = 0;
        for (auto stratum : all_strata) {
            for (int j = 0; j < n_windows; j++) {

                xt::col(all_windows, j) += xt::sum(xt::view(stratum, xt::all(),
                                                            xt::range(j * window_size,
                                                                      (j + 1) *
                                                                      window_size -
                                                                      stra_i)),
                                                   1);
            }
            stra_i++;
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
//    if (print_time) {
//        cout << "Time 1:" << duration1.count() << endl
//             << "Time 2:" << duration2.count() << endl;
//    }
    cout << "Time 1:" << duration1.count() << endl
         << "Time 2:" << duration2.count() << endl;
    ofstream fout;
    fout.open("time_test.txt");
    fout<<"Time 1: "<< duration1.count() << endl<< "Time 2:" << duration2.count() <<
    endl;
    fout.close();
    double tout = duration1.count();
    return make_pair(distance_mat,tout) ;
}


