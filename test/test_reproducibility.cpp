////
//// Created by Bean Juice on 06/01/2021.
////
//
//#include "xtensor/xio.hpp"
//#include "../load/ContactMaps.h"
//#include "../embedding/reproducibility.h"
//
//using namespace std;
//
////}
//double innerP(vector<xt::xarray<double>> all_strata) {
//    xt::xarray<double> pair_dis;
//    double time1;
//    tie(pair_dis, time1) =
//            pairwise_distance(all_strata, "inner_product");
//    cout << "pairwise dis: " << pair_dis << endl;
//    return time1;
//}
//
//double selfishP(vector<xt::xarray<double>> all_strata) {
//    xt::xarray<double> pair_dis;
//    double time1;
//    tie(pair_dis, time1) =
//            pairwise_distance(all_strata, "selfish");
//    cout << "pairwise dis: " << pair_dis << endl;
//    return time1;
//}
//int main() {
//    vector<string> fileLst{"../test/data/cell_03", "../test/data/cell_01",
//                           "../test/data/cell_02"};
//    vector<string> operation{"convolution"};
//    scHiCs y = scHiCs(fileLst, "mm9", 500000, 3, 4000000, true, "except Y",
//                      "shortest_score",
//                      10, true,
//                      operation);
//    vector<xt::xarray<double>> all_strata = y.get_strata()["chr8"];
//    double tsum = innerP(all_strata);
//    //xt::xarray<double> pair_dis;double time1; tie(pair_dis,time1)= pairwise_distance
//    //      (all_strata, "hicrep");
//    cout << "time: " << tsum;
//}
//
