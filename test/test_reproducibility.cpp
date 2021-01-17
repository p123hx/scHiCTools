//
// Created by Bean Juice on 06/01/2021.
//

#include "xtensor/xio.hpp"
#include "../load/ContactMaps.h"
#include "../embedding/reproducibility.h"

using namespace std;

int main() {
    vector<string> fileLst{"../test/data/cell_03", "../test/data/cell_01",
                           "../test/data/cell_02"};
    vector<string> operation{"convolution"};
    scHiCs y = scHiCs(fileLst, "mm9", 100000, 3, 4000000, false, "except Y",
                      "shortest_score",
                      10, true,
                      operation);
    vector<xt::xarray<double>> all_strata = y.get_strata()["chr8"];
    xt::xarray<double> pair_dis = pairwise_distance(all_strata, "hicrep");

    cout << "pairwise dis: " << pair_dis << endl;

}