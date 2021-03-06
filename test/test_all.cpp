#include "xtensor/xio.hpp"
#include "../load/ContactMaps.h"
#include "../embedding/reproducibility.h"
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>

using namespace std;

vector<double> innerP(vector<xt::xarray<double>> all_strata) {
    xt::xarray<double> pair_dis;
    vector<double> times;
    tie(pair_dis, times) =
            pairwise_distance(all_strata, "inner_product");
//    cout << "pairwise dis: " << pair_dis << endl;
    return times;
}

vector<double> fastHicP(vector<xt::xarray<double>> all_strata) {
    xt::xarray<double> pair_dis;
    vector<double> times;
    tie(pair_dis, times) =
            pairwise_distance(all_strata, "hicrep");
//   cout << "pairwise dis: " << pair_dis << endl;
    return times;
}

double selfishP(vector<xt::xarray<double>> all_strata) {
    xt::xarray<double> pair_dis;
    vector<double> times;
    tie(pair_dis, times) =
            pairwise_distance(all_strata, "selfish");
//    cout << "pairwise dis: " << pair_dis << endl;
    return times[0];
}

double oldHicP(vector<xt::xarray<double>> all_strata) {
    xt::xarray<double> pair_dis;
    vector<double> times;
    tie(pair_dis, times) =
            pairwise_distance(all_strata, "old_hicrep");
//    cout << "pairwise dis: " << pair_dis << endl;
    return times[0];
}

vector<string> f1000() {
    return vector<string>{"../../Nagano/1CDX_cells/1CDX1.1/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.185/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.281/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.38/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.46/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.117/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.202/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.294/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.377/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.465/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.108/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.182/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.263/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.352/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.68/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.154/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.237/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.312/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.392/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.468/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.101/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.186/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.283/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.381/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.464/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.12/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.203/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.295/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.382/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.466/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.11/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.183/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.264/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.353/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.72/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.155/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.24/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.313/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.393/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.47/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.102/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.187/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.284/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.383/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.466/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.121/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.204/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.296/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.383/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.467/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.111/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.185/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.265/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.354/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.73/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.156/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.241/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.314/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.394/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.472/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.103/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.191/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.285/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.384/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.468/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.122/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.205/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.297/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.384/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.468/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.112/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.186/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.266/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.355/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.74/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.157/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.242/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.315/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.396/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.473/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.104/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.192/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.286/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.385/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.47/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.123/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.206/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.3/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.385/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.47/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.113/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.187/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.267/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.356/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.75/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.158/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.243/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.316/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.397/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.474/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.105/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.193/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.287/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.386/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.472/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.124/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.207/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.301/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.387/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.472/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.114/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.188/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.27/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.357/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.76/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.16/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.244/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.317/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.4/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.475/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.106/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.194/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.293/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.387/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.474/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.125/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.21/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.302/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.391/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.473/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.115/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.191/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.271/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.358/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.81/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.161/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.245/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.32/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.403/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.476/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.107/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.195/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.294/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.388/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.475/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.126/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.211/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.303/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.392/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.474/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.117/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.192/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.272/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.36/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.82/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.162/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.246/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.321/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.404/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.477/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.108/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.196/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.295/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.391/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.476/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.13/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.212/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.304/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.393/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.475/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.12/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.193/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.273/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.362/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.83/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.163/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.247/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.323/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.405/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.48/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.111/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.197/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.296/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.392/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.477/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.131/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.213/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.305/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.394/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.476/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.121/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.194/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.274/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.363/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.85/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.164/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.248/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.324/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.406/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.482/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.113/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.2/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.297/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.393/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.482/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.132/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.214/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.306/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.395/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.477/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.122/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.195/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.275/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.364/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.86/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.165/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.25/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.325/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.407/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.483/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.114/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.201/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.3/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.394/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.483/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.133/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.215/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.307/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.396/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.48/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.123/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.196/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.276/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.365/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.87/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.166/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.251/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.326/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.41/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.484/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.115/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.202/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.303/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.395/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.484/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.134/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.216/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.31/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.4/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.481/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.124/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.197/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.277/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.366/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.91/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.167/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.252/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.327/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.411/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.485/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.116/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.203/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.304/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.396/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.487/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.135/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.217/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.312/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.401/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.482/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.125/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.2/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.282/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.368/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.92/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.17/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.253/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.33/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.412/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.486/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.117/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.204/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.305/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.397/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.5/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.136/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.22/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.313/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.402/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.483/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.126/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.201/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.283/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.372/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.93/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.171/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.254/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.332/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.413/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.487/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.12/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.205/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.306/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.398/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.51/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.137/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.221/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.314/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.403/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.486/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.127/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.202/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.284/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.373/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.94/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.172/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.255/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.333/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.414/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.5/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.121/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.207/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.307/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.4/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.52/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.14/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.222/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.315/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.404/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.487/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.128/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.203/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.285/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.382/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.95/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.173/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.256/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.334/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.415/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.51/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.122/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.211/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.312/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.402/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.53/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.141/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.224/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.316/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.405/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.488/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.13/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.204/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.287/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.383/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.96/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.174/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.257/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.335/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.416/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.52/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.123/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.213/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.313/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.403/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.54/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.142/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.225/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.317/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.407/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.5/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.131/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.205/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.293/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.384/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.97/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.175/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.258/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.336/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.417/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.53/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.124/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.214/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.314/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.404/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.56/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.144/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.227/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.322/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.408/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.51/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.132/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.207/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.294/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.387/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.1/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.176/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.26/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.337/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.418/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.54/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.125/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.215/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.315/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.405/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.61/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.145/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.23/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.323/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.41/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.52/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.133/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.21/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.295/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.392/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.101/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.177/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.261/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.34/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.42/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.55/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.126/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.216/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.317/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.406/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.62/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.146/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.232/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.324/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.412/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.53/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.134/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.211/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.296/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.395/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.102/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.178/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.262/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.341/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.421/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.57/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.127/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.22/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.32/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.407/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.63/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.147/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.233/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.325/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.413/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.54/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.135/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.212/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.297/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.406/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.103/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.181/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.263/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.342/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.422/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.6/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.128/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.221/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.322/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.41/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.64/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.15/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.234/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.326/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.414/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.55/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.136/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.213/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.3/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.411/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.104/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.182/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.264/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.343/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.423/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.61/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.13/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.222/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.323/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.411/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.66/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.151/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.235/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.327/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.415/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.56/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.137/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.214/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.301/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.412/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.105/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.183/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.265/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.344/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.424/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.62/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.132/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.223/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.324/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.412/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.67/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.154/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.236/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.33/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.416/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.57/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.138/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.215/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.302/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.42/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.106/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.185/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.266/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.345/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.425/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.63/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.134/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.224/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.325/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.413/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.71/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.155/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.237/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.331/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.417/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.58/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.14/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.216/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.303/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.424/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.107/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.186/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.267/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.346/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.426/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.64/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.136/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.226/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.326/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.414/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.72/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.156/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.24/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.332/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.42/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.6/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.141/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.217/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.304/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.426/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.11/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.187/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.268/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.347/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.428/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.65/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.138/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.23/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.327/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.415/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.73/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.157/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.241/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.333/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.421/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.62/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.142/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.218/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.305/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.428/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.111/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.191/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.27/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.348/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.43/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.66/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.14/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.232/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.33/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.416/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.74/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.16/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.242/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.334/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.422/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.63/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.143/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.22/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.306/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.43/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.113/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.192/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.272/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.35/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.431/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.67/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.141/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.233/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.333/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.417/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.75/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.162/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.243/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.335/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.423/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.64/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.144/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.221/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.307/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.437/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.114/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.193/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.273/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.352/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.432/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.71/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.142/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.234/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.334/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.42/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.76/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.163/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.244/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.336/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.424/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.65/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.145/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.222/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.31/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.438/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.115/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.194/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.274/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.353/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.433/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.72/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.143/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.235/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.335/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.421/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.77/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.164/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.245/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.337/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.425/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.67/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.146/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.223/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.311/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.443/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.116/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.195/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.275/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.354/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.434/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.73/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.144/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.236/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.337/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.422/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.81/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.165/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.25/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.34/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.426/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.72/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.147/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.224/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.312/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.447/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.117/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.197/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.276/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.355/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.435/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.74/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.145/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.24/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.34/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.423/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.82/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.166/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.254/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.342/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.427/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.73/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.15/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.226/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.313/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.452/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.12/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.2/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.277/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.356/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.436/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.75/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.146/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.242/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.342/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.424/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.83/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.167/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.255/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.343/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.43/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.74/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.151/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.227/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.314/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.455/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.121/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.201/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.278/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.357/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.437/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.76/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.147/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.244/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.343/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.425/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.85/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.17/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.256/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.344/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.432/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.75/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.152/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.23/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.315/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.46/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.124/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.202/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.281/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.36/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.438/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.77/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.151/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.245/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.344/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.426/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.86/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.171/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.257/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.345/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.433/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.76/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.153/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.232/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.316/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.466/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.125/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.203/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.282/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.363/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.441/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.81/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.152/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.25/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.345/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.427/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.87/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.172/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.26/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.346/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.435/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.77/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.154/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.233/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.317/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.47/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.126/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.204/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.283/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.364/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.442/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.82/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.153/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.252/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.346/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.43/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.92/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.173/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.262/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.348/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.436/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.82/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.155/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.234/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.32/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.473/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.127/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.205/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.284/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.365/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.443/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.83/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.154/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.253/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.347/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.432/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.93/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.174/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.263/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.35/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.437/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.83/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.156/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.235/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.322/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.475/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.13/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.206/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.285/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.366/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.444/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.84/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.155/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.254/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.35/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.433/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.94/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.175/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.264/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.352/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.44/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.84/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.157/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.236/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.323/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.483/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.132/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.207/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.286/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.367/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.446/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.85/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.156/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.255/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.352/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.434/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.95/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.176/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.265/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.353/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.441/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.85/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.158/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.237/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.324/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.484/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.133/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.212/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.287/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.37/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.447/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.86/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.157/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.257/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.353/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.435/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.96/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.177/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.266/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.354/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.442/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.86/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.16/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.24/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.326/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.486/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.134/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.214/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.291/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.371/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.448/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.87/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.161/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.258/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.354/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.436/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.1/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.181/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.267/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.355/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.443/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.87/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.161/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.242/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.327/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.5/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.135/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.216/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.292/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.372/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.45/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.91/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.162/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.26/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.355/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.437/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.101/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.182/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.27/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.356/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.444/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.91/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.163/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.243/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.33/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.51/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.136/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.217/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.293/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.373/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.452/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.92/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.164/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.262/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.357/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.44/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.102/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.183/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.272/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.357/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.445/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.92/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.164/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.244/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.332/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.52/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.137/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.22/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.294/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.374/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.453/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.93/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.165/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.263/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.36/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.442/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.103/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.184/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.273/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.36/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.446/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.93/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.165/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.245/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.333/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.53/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.14/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.221/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.295/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.375/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.454/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.94/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.166/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.264/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.362/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.443/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.104/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.185/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.274/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.362/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.447/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.94/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.166/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.246/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.334/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.54/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.142/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.222/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.296/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.376/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.455/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.95/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.167/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.265/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.363/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.444/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.105/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.186/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.275/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.363/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.45/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.95/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.168/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.248/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.335/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.55/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.143/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.223/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.297/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.377/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.456/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.97/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.17/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.266/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.364/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.445/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.106/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.187/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.276/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.364/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.452/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.96/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.17/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.25/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.337/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.56/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.144/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.224/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.3/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.38/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.457/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.171/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.267/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.365/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.447/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.107/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.191/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.277/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.365/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.454/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.97/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.171/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.251/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.34/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.57/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.145/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.225/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.301/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.381/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.458/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.172/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.27/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.366/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.45/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.108/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.192/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.282/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.366/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.455/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.1/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.172/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.252/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.342/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.6/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.146/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.226/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.302/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.382/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.46/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.174/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.272/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.37/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.452/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.11/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.193/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.283/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.367/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.456/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.102/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.173/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.253/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.343/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.61/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.147/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.23/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.303/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.383/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.462/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.175/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.273/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.371/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.453/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.112/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.194/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.284/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.37/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.457/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.103/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.174/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.254/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.344/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.62/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.148/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.232/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.304/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.384/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.463/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.176/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.274/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.372/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.454/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.113/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.195/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.285/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.372/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.46/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.104/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.175/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.255/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.345/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.64/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.15/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.233/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.305/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.385/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.464/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.177/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.275/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.373/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.455/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.114/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.196/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.286/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.374/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.461/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.105/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.176/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.256/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.346/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.65/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.151/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.234/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.307/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.386/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.465/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.182/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.276/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.374/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.456/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.115/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.197/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.287/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.375/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.463/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.106/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.177/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.257/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.347/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.66/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.152/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.235/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.31/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.387/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.466/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.183/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.277/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.375/new_adj",
                          "../../Nagano/1CDX_cells/1CDX1.457/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.116/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.2/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.293/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.376/new_adj",
                          "../../Nagano/1CDX_cells/1CDX2.464/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.107/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.181/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.262/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.35/new_adj",
                          "../../Nagano/1CDX_cells/1CDX3.67/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.153/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.236/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.311/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.391/new_adj",
                          "../../Nagano/1CDX_cells/1CDX4.467/new_adj"};
}

/*void toolN(int n){
    vector<string> fileLst1000 = f1000();
    vector<string>fileLstN (fileLst1000.begin(),fileLst1000.begin()+n);
    vector<string> operation{"convolution"};
    scHiCs y = scHiCs(fileLstN, "mm9", 500000, 3, 4000000, true, "except Y",
                         "shortest_score",
                         10, false,
                         operation);
    vector<string> chrs{"chr1", "chr2", "chrX", "chr3", "chr4", "chr5", "chr6", "chr7",
                        "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                        "chr15", "chr16", "chr17", "chr18", "chr19"}; //since "except Y"
    double inner100 = 0.0, selfish100 = 0.0, inner100t1 = 0.0, inner100t2 = 0.0,
            fast100 = 0.0;
    for(int i=10;i>0;i--){
        for (string s:chrs) {
            cout << "\n" << s << ":\n";
            vector<xt::xarray<double>> chr = y.get_strata()[s];
            fast100 += (fastHicP(chr) / 100000);
            vector<double> tmpD = innerP(chr);
            inner100 += (tmpD[0] / 100000);
            inner100t1 += (tmpD[1] / 100000);
            inner100t2 += (tmpD[2] / 100000);
            selfish100 += (selfishP(chr) / 100000);
        }
    }

    cout << "inner 100 cells:\n t1: " << inner100t1 << " t2: " << inner100t2
         << " total: " << inner100 << " in seconds\n";
    cout << "time1 + time2 selfish 100 cells: " << selfish100 << " in "
                                                                 "seconds\n"
         << "fast 100 new:"<<fast100<<endl;                                                  "milliseconds\n";
    string outF = to_string(n)+"out.txt";
    ofstream fout(outF);
    fout << "inner 100 cells:\n t1: " << inner100t1 << " t2: " << inner100t2
         << " total: " << inner100 << " in "

                                      "seconds\n";
    fout << "time1 + time2 selfish 100 cells: " << selfish100 << " in "
                                                                 "seconds\n"
         << "fast 100 new:"<<fast100<<endl;
    fout.close();
}
void test(){
    vector<string> fileLst{"../NaganoPartial/1CDX1.1/new_adj","../NaganoPartial/1CDX1"
                                                              ".2/new_adj"};
    vector<string> operation{"convolution"};
    scHiCs y = scHiCs(fileLst, "mm9", 500000, 3, 4000000, true, "except Y",
                      "shortest_score",
                      10, false,
                      operation);
    vector<string> chrs{"chr1", "chr2", "chrX", "chr3", "chr4", "chr5", "chr6", "chr7",
                        "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                        "chr15", "chr16", "chr17", "chr18", "chr19"}; //since "except Y"
    double tsum = 0.0;
    for(int i=10;i>0;i--){
        for (string s:chrs) {
            cout << "\n" << s << ":\n";
            vector<xt::xarray<double>> chr = y.get_strata()[s];
            tsum += (selfishP(chr) / 100000);
        }
    }

    cout << "Total average fastHiCrep test " << tsum<< " in milliseconds\n";
}*/
void allNew(int n) {
    vector<string> fileLst1000 = f1000();
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    shuffle(fileLst1000.begin(), fileLst1000.end(), std::default_random_engine(seed));

    vector<string> fileLstN(fileLst1000.begin(), fileLst1000.begin() + n);
    vector<string> operation{"convolution"};
    ofstream allF("allF.csv", std::ios::app);
    ofstream allI("allI.csv", std::ios::app);
    scHiCs y = scHiCs(fileLstN, "mm9", 500000, 3, 4000000, true, "except Y",
                      "shortest_score",
                      10, true,
                      operation);
    vector<string> chrs{"chr1", "chr2", "chrX", "chr3", "chr4", "chr5", "chr6", "chr7",
                        "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                        "chr15", "chr16", "chr17", "chr18", "chr19"}; //since "except Y"
    double innerT = 0.0, innerTt1 = 0.0, innerTt2 = 0.0, innerTt3 = 0.0, innerTt4 = 0.0,
            fastT = 0.0, fastTt1 = 0.0, fastTt2 = 0.0, fastTt3 = 0.0, fastTt4 = 0.0,
            innerTM = 0.0, innerTt1M = 0.0, innerTt2M = 0.0, innerTt3M = 0.0, innerTt4M = 0.0,
            fastTM = 0.0, fastTt1M = 0.0, fastTt2M = 0.0, fastTt3M = 0.0, fastTt4M = 0.0;
    string outF = to_string(n) + "out.txt";
    ofstream fout(outF);
    for (int i = 11; i > 0; i--) {
        double tinnerT = 0.0, tinnerTt1 = 0.0, tinnerTt2 = 0.0, tinnerTt3 = 0.0,
                tinnerTt4 = 0.0,
                tfastT = 0.0, tfastTt1 = 0.0, tfastTt2 = 0.0, tfastTt3 = 0.0,
                tfastTt4 = 0.0;
        cout << "Set: " << i << endl;
        for (string s:chrs) {
            vector<xt::xarray<double>> chr = y.get_strata()[s];
            vector<double> tmpD = innerP(chr);
            vector<double> tmpI = fastHicP(chr);
            tinnerT += (tmpI[0]) / 1000.0;
            tinnerTt1 += (tmpI[1]) / 1000.0;
            tinnerTt2 += (tmpI[2]) / 1000.0;
            tinnerTt3 += (tmpI[3]) / 1000.0;
            tinnerTt4 += (tmpI[4]) / 1000.0;
            tfastT += (tmpD[0]) / 1000.0;
            tfastTt1 += (tmpD[1]) / 1000.0;
            tfastTt2 += (tmpD[2]) / 1000.0;
            tfastTt3 += (tmpD[3]) / 1000.0;
            tfastTt4 += (tmpD[4]) / 1000.0;
        }
        innerTM = max(innerTM, tinnerT);
        innerTt1M = max(innerTt1M, tinnerTt1);
        innerTt2M = max(innerTt2M, tinnerTt2);
        innerTt3M = max(innerTt3M, tinnerTt3);
        innerTt4M = max(innerTt4M, tinnerTt4);
        fastTM = max(fastTM, tfastT);
        fastTt1M = max(fastTt1M, tfastTt1);
        fastTt2M = max(fastTt2M, tfastTt2);
        fastTt3M = max(fastTt3M, tfastTt3);
        fastTt4M = max(fastTt4M, tfastTt4);
        innerT += tinnerT;
        innerTt1 += tinnerTt1;
        innerTt2 += tinnerTt2;
        innerTt3 += tinnerTt3;
        innerTt4 += tinnerTt4;
        fastT += tfastT;
        fastTt1 += tfastTt1;
        fastTt2 += tfastTt2;
        fastTt3 += tfastTt3;
        fastTt4 += tfastTt4;
        cout << n << ": fast set: " << tfastT << " t1: " << tfastTt1 << " t2: " <<
             tfastTt2 << " t3: " << tfastTt3 << " t4: " << tfastTt4 << endl
             << " inner set: " << tinnerT << " t1: " << tinnerTt1 << " t2: " <<
             tinnerTt2 << " t3: " << tinnerTt3 << " t4: " << tinnerTt4 << endl;
        fout << n << ": fast set: " << tfastT << "; t1: " << tinnerTt1 << " t2: "
             << tinnerTt2
             << " set_total: " << tinnerT << endl;
    }
    innerT -= innerTM;
    innerTt1 -= innerTt1M;
    innerTt2 -= innerTt2M;
    innerTt3 -= innerTt3M;
    innerTt4 -= innerTt4M;
    fastT -= fastTM;
    fastTt1 -= fastTt1M;
    fastTt2 -= fastTt2M;
    fastTt3 -= fastTt3M;
    fastTt4 -= fastTt4M;

    innerT /= 10.0;
    innerTt1 /= 10.0;
    innerTt2 /= 10.0;
    innerTt3 /= 10.0;
    innerTt4 /= 10.0;
    fastT /= 10.0;
    fastTt1 /= 10.0;
    fastTt2 /= 10.0;
    fastTt3 /= 10.0;
    fastTt4 /= 10.0;
    cout << "fast total:" << fastT << " t1: " << fastTt1 << " t2: "
         << fastTt2 << " t3: " << fastTt3 << " t4: " << fastTt4 << endl
         << "inner total: " << innerT << " t1: " << innerTt1 << " t2: "
         << innerTt2 << " t3: " << innerTt3 << " t4: " << innerTt4
         << endl;


    fout << "fast total:" << fastT << " t1: " << fastTt1 << " t2: "
         << fastTt2 << " t3: " << fastTt3 << " t4: " << fastTt4 << endl
         << "inner total: " << innerT << " t1: " << innerTt1 << " t2: "
         << innerTt2 << " t3: " << innerTt3 << " t4: " << innerTt4
         << endl;
    fout.close();
    allF << n << "," << fastT << "," << fastTt1 << "," << fastTt2 << "," << fastTt3 << "," << fastTt4 << endl;
    allI << n << "," << innerT << "," << innerTt1 << "," << innerTt2 << "," << innerTt3 << ","
                                                                                           "" << innerTt4 << endl;
    for (int cellC = n; cellC < 1000; cellC += n) {
        vector<string> fileLst100(fileLst1000.begin() + cellC, fileLst1000.begin() +
                                                               n + cellC);
        y.load100(fileLst100, "mm9", 500000, 3, 4000000, true, "except Y",
                  "shortest_score",
                  10, true,
                  operation);
        innerT = 0.0, innerTt1 = 0.0, innerTt2 = 0.0, innerTt3 = 0.0, innerTt4 = 0.0,
        fastT = 0.0, fastTt1 = 0.0, fastTt2 = 0.0, fastTt3 = 0.0, fastTt4 = 0.0,
        innerTM = 0.0, innerTt1M = 0.0, innerTt2M = 0.0, innerTt3M = 0.0, innerTt4M = 0.0,
        fastTM = 0.0, fastTt1M = 0.0, fastTt2M = 0.0, fastTt3M = 0.0, fastTt4M = 0.0;
        outF = to_string(cellC + n) + "out.txt";
        ofstream fout(outF);
        for (int i = 11; i > 0; i--) {
            double tinnerT = 0.0, tinnerTt1 = 0.0, tinnerTt2 = 0.0, tinnerTt3 = 0.0,
                    tinnerTt4 = 0.0,
                    tfastT = 0.0, tfastTt1 = 0.0, tfastTt2 = 0.0, tfastTt3 = 0.0,
                    tfastTt4 = 0.0;
            cout << "Set: " << i << endl;
            for (string s:chrs) {
                vector<xt::xarray<double>> chr = y.get_strata()[s];
                vector<double> tmpD = innerP(chr);
                vector<double> tmpI = fastHicP(chr);
                tinnerT += (tmpI[0]) / 1000.0;
                tinnerTt1 += (tmpI[1]) / 1000.0;
                tinnerTt2 += (tmpI[2]) / 1000.0;
                tinnerTt3 += (tmpI[3]) / 1000.0;
                tinnerTt4 += (tmpI[4]) / 1000.0;
                tfastT += (tmpD[0]) / 1000.0;
                tfastTt1 += (tmpD[1]) / 1000.0;
                tfastTt2 += (tmpD[2]) / 1000.0;
                tfastTt3 += (tmpD[3]) / 1000.0;
                tfastTt4 += (tmpD[4]) / 1000.0;
            }
            innerTM = max(innerTM, tinnerT);
            innerTt1M = max(innerTt1M, tinnerTt1);
            innerTt2M = max(innerTt2M, tinnerTt2);
            innerTt3M = max(innerTt3M, tinnerTt3);
            innerTt4M = max(innerTt4M, tinnerTt4);
            fastTM = max(fastTM, tfastT);
            fastTt1M = max(fastTt1M, tfastTt1);
            fastTt2M = max(fastTt2M, tfastTt2);
            fastTt3M = max(fastTt3M, tfastTt3);
            fastTt4M = max(fastTt4M, tfastTt4);
            innerT += tinnerT;
            innerTt1 += tinnerTt1;
            innerTt2 += tinnerTt2;
            innerTt3 += tinnerTt3;
            innerTt4 += tinnerTt4;
            fastT += tfastT;
            fastTt1 += tfastTt1;
            fastTt2 += tfastTt2;
            fastTt3 += tfastTt3;
            fastTt4 += tfastTt4;
            cout << cellC + n << ": fast set: " << tfastT << " t1: " << tfastTt1 << " t2: " <<
                 tfastTt2 << " t3: " << tfastTt3 << " t4: " << tfastTt4 << endl
                 << " inner set: " << tinnerT << " t1: " << tinnerTt1 << " t2: " <<
                 tinnerTt2 << " t3: " << tinnerTt3 << " t4: " << tinnerTt4 << endl;
            fout << cellC + n << ": fast set: " << tfastT << "; t1: " << tinnerTt1 << " t2: "
                 << tinnerTt2
                 << " set_total: " << tinnerT << endl;
        }
        innerT -= innerTM;
        innerTt1 -= innerTt1M;
        innerTt2 -= innerTt2M;
        innerTt3 -= innerTt3M;
        innerTt4 -= innerTt4M;
        fastT -= fastTM;
        fastTt1 -= fastTt1M;
        fastTt2 -= fastTt2M;
        fastTt3 -= fastTt3M;
        fastTt4 -= fastTt4M;

        innerT /= 10.0;
        innerTt1 /= 10.0;
        innerTt2 /= 10.0;
        innerTt3 /= 10.0;
        innerTt4 /= 10.0;
        fastT /= 10.0;
        fastTt1 /= 10.0;
        fastTt2 /= 10.0;
        fastTt3 /= 10.0;
        fastTt4 /= 10.0;
        cout << "fast total:" << fastT << " t1: " << fastTt1 << " t2: "
             << fastTt2 << " t3: " << fastTt3 << " t4: " << fastTt4 << endl
             << "inner total: " << innerT << " t1: " << innerTt1 << " t2: "
             << innerTt2 << " t3: " << innerTt3 << " t4: " << innerTt4
             << endl;


        fout << "fast total:" << fastT << " t1: " << fastTt1 << " t2: "
             << fastTt2 << " t3: " << fastTt3 << " t4: " << fastTt4 << endl
             << "inner total: " << innerT << " t1: " << innerTt1 << " t2: "
             << innerTt2 << " t3: " << innerTt3 << " t4: " << innerTt4
             << endl;
        fout.close();
        allF << cellC + n << "," << fastT << "," << fastTt1 << "," << fastTt2 << "," << fastTt3 << ","
                                                                                                   "" << fastTt4 << endl;
        allI << cellC + n << "," << innerT << "," << innerTt1 << "," << innerTt2 << "," << innerTt3 << ","
                                                                                                       "" << innerTt4 << endl;
    }
    allF.close();
    allI.close();
}

void testNew() {
    vector<string> fileLst{"../test/data/cell_03"};
    vector<string> operation{"convolution"};

    scHiCs y = scHiCs(fileLst, "mm9", 500000, 3, 4000000, true, "except Y",
                      "shortest_score",
                      10, true,
                      operation);
    vector<string> file1{"../test/data/cell_01"};
    vector<string> file2{"../test/data/cell_02"};
    y.load100(file1, "mm9", 500000, 3, 4000000, true, "except Y",
              "shortest_score",
              10, true,
              operation);
    y.load100(file2, "mm9", 500000, 3, 4000000, true, "except Y",
              "shortest_score",
              10, true,
              operation);
    vector<xt::xarray<double>> chr = y.get_strata()["chr1"];
    vector<double> tmpD = innerP(chr);
}

void toolNew(int n) {
    vector<string> fileLst1000 = f1000();
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    shuffle(fileLst1000.begin(), fileLst1000.end(), std::default_random_engine(seed));

    vector<string> fileLstN(fileLst1000.begin(), fileLst1000.begin() + n);
    vector<string> operation{"convolution"};
    ofstream allF("allF.csv", std::ios::app);
    ofstream allI("allI.csv", std::ios::app);
    scHiCs y = scHiCs(fileLstN, "mm9", 500000, 3, 4000000, true, "except Y",
                      "shortest_score",
                      10, true,
                      operation);
    vector<string> chrs{"chr1", "chr2", "chrX", "chr3", "chr4", "chr5", "chr6", "chr7",
                        "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                        "chr15", "chr16", "chr17", "chr18", "chr19"}; //since "except Y"
    double innerT = 0.0, innerTt1 = 0.0, innerTt2 = 0.0, innerTt3 = 0.0, innerTt4 = 0.0,
            fastT = 0.0, fastTt1 = 0.0, fastTt2 = 0.0, fastTt3 = 0.0, fastTt4 = 0.0,
            innerTM = 0.0, innerTt1M = 0.0, innerTt2M = 0.0, innerTt3M = 0.0, innerTt4M = 0.0,
            fastTM = 0.0, fastTt1M = 0.0, fastTt2M = 0.0, fastTt3M = 0.0, fastTt4M = 0.0;
    string outF = to_string(n) + "out.txt";
    ofstream fout(outF);
    for (int i = 11; i > 0; i--) {
        double tinnerT = 0.0, tinnerTt1 = 0.0, tinnerTt2 = 0.0, tinnerTt3 = 0.0,
                tinnerTt4 = 0.0,
                tfastT = 0.0, tfastTt1 = 0.0, tfastTt2 = 0.0, tfastTt3 = 0.0,
                tfastTt4 = 0.0;
        cout << "Set: " << i << endl;
        for (string s:chrs) {
            vector<xt::xarray<double>> chr = y.get_strata()[s];
            vector<double> tmpD = innerP(chr);
            vector<double> tmpI = fastHicP(chr);
            tinnerT += (tmpI[0]) / 1000.0;
            tinnerTt1 += (tmpI[1]) / 1000.0;
            tinnerTt2 += (tmpI[2]) / 1000.0;
            tinnerTt3 += (tmpI[3]) / 1000.0;
            tinnerTt4 += (tmpI[4]) / 1000.0;
            tfastT += (tmpD[0]) / 1000.0;
            tfastTt1 += (tmpD[1]) / 1000.0;
            tfastTt2 += (tmpD[2]) / 1000.0;
            tfastTt3 += (tmpD[3]) / 1000.0;
            tfastTt4 += (tmpD[4]) / 1000.0;
        }
        innerTM = max(innerTM, tinnerT);
        innerTt1M = max(innerTt1M, tinnerTt1);
        innerTt2M = max(innerTt2M, tinnerTt2);
        innerTt3M = max(innerTt3M, tinnerTt3);
        innerTt4M = max(innerTt4M, tinnerTt4);
        fastTM = max(fastTM, tfastT);
        fastTt1M = max(fastTt1M, tfastTt1);
        fastTt2M = max(fastTt2M, tfastTt2);
        fastTt3M = max(fastTt3M, tfastTt3);
        fastTt4M = max(fastTt4M, tfastTt4);
        innerT += tinnerT;
        innerTt1 += tinnerTt1;
        innerTt2 += tinnerTt2;
        innerTt3 += tinnerTt3;
        innerTt4 += tinnerTt4;
        fastT += tfastT;
        fastTt1 += tfastTt1;
        fastTt2 += tfastTt2;
        fastTt3 += tfastTt3;
        fastTt4 += tfastTt4;
        cout << n << ": fast set: " << tfastT << " t1: " << tfastTt1 << " t2: " <<
             tfastTt2 << " t3: " << tfastTt3 << " t4: " << tfastTt4 << endl
             << " inner set: " << tinnerT << " t1: " << tinnerTt1 << " t2: " <<
             tinnerTt2 << " t3: " << tinnerTt3 << " t4: " << tinnerTt4 << endl;
        fout << n << ": fast set: " << tfastT << "; t1: " << tinnerTt1 << " t2: "
             << tinnerTt2
             << " set_total: " << tinnerT << endl;
    }
    innerT -= innerTM;
    innerTt1 -= innerTt1M;
    innerTt2 -= innerTt2M;
    innerTt3 -= innerTt3M;
    innerTt4 -= innerTt4M;
    fastT -= fastTM;
    fastTt1 -= fastTt1M;
    fastTt2 -= fastTt2M;
    fastTt3 -= fastTt3M;
    fastTt4 -= fastTt4M;

    innerT /= 10.0;
    innerTt1 /= 10.0;
    innerTt2 /= 10.0;
    innerTt3 /= 10.0;
    innerTt4 /= 10.0;
    fastT /= 10.0;
    fastTt1 /= 10.0;
    fastTt2 /= 10.0;
    fastTt3 /= 10.0;
    fastTt4 /= 10.0;
    cout << "fast total:" << fastT << " t1: " << fastTt1 << " t2: "
         << fastTt2 << " t3: " << fastTt3 << " t4: " << fastTt4 << endl
         << "inner total: " << innerT << " t1: " << innerTt1 << " t2: "
         << innerTt2 << " t3: " << innerTt3 << " t4: " << innerTt4
         << endl;


    fout << "fast total:" << fastT << " t1: " << fastTt1 << " t2: "
         << fastTt2 << " t3: " << fastTt3 << " t4: " << fastTt4 << endl
         << "inner total: " << innerT << " t1: " << innerTt1 << " t2: "
         << innerTt2 << " t3: " << innerTt3 << " t4: " << innerTt4
         << endl;
    fout.close();
    allF << n << "," << fastT << "," << fastTt1 << "," << fastTt2 << "," << fastTt3 << "," << fastTt4 << endl;
    allI << n << "," << innerT << "," << innerTt1 << "," << innerTt2 << "," << innerTt3 << ","
                                                                                           "" << innerTt4 << endl;
    allF.close();
    allI.close();
}

void test() {
    vector<string> fileLst{"../NaganoPartial/1CDX1.1/new_adj", "../NaganoPartial/1CDX1"
                                                               ".2/new_adj", "../NaganoPartial/1CDX1.3/new_adj"};
    vector<string> operation{"convolution"};

    scHiCs y = scHiCs(fileLst, "mm9", 500000, 3, 4000000, true, "except Y",
                      "shortest_score",
                      10, true,
                      operation);
    vector<xt::xarray<double>> chr = y.get_strata()["chr1"];
    vector<double> tmpD = innerP(chr);
}

int main() {
    vector<int> v{100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
    for (int i:v) {
        allNew(i);
    }
}


