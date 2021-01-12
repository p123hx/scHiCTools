//
// Created by Bean Juice on 06/01/2021.
//

#include <utility>
#include "xtensor/xio.hpp"
#include "xtensor/xview.hpp"
#include "xtensor-blas/xlinalg.hpp"
#include "../load/ContactMaps.h"
using namespace std;
int main(){
vector<string> fileLst {"data/cell_03","data/cell_01","data/cell_02"};
vector<string> operation {"convolution"};
scHiCs y = scHiCs(fileLst,"mm9",100000,true,false, "all","shortest_score", 10,false,
                  operation,0,0,0,false,false,0);
}