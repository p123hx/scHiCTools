//
// Created by Bean Juice on 16/12/2020.
//

#ifndef SCHICTOOLS_LOAD_HIC_FILE_H
#define SCHICTOOLS_LOAD_HIC_FILE_H

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

//pair<set<string>, map<string, int>>
//get_chromosome_lengths(const string &ref_str, const string &chromosomes,
//                       set<string> chrom_set, int res);
pair<set<string>, map<string, int>>
get_chromosome_lengths(const string &ref_str, const string &chromosomes,
                       int res);

//get_chromosome_lengths(map<string, int> ref_map, string chromosomes,
//                       set<string> chrom_set, int res);
pair<set<string>, map<string, int>>
get_chromosome_lengths(map<string, int> ref_map, string chromosomes,
                       int res);

pair<xt::xarray<double>, vector<xt::xarray<double>>> load_HiC(string file, map<string,
        int> genome_length,
                                                              string
                                                              format = "",
                                                              int custom_format =
                                                              0,
                                                              int header = 0,
                                                              string chromosome = "",
                                                              int resolution = 10000,
                                                              bool resolution_adjust = true,
                                                              double map_filter = 0.,
                                                              bool sparse = false,
                                                              bool gzip = false,
                                                              int keep_n_strata = 0,
                                                              vector<string> operations = vector<string>());

#endif //SCHICTOOLS_LOAD_HIC_FILE_H

