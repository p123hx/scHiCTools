//
// Created by Bean Juice on 06/01/2021.
//

#include <utility>
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include "xtensor-blas/xlinalg.hpp"
using namespace std;
int main(){
xt::xarray <double > mat;
vector<xt::xarray <double>> strata;

tie(mat,strata) =      scHiCs(['data/cell_03','data/cell_01','data/cell_02'],
reference_genome='mm9', resolution=100000,
max_distance=4000000, format='shortest_score',
adjust_resolution=False, chromosomes='except Y',
operations=['convolution'], kernel_shape=3, keep_n_strata=10,
                            store_full_map=True)
}