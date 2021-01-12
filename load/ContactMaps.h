//
// Created by Bean Juice on 16/12/2020.
//

#ifndef SCHICTOOLS_CONTACTMAPS_H
#define SCHICTOOLS_CONTACTMAPS_H

#include <string>
#include "ContactMaps.h"
#include "load_hic_file.h"
#include <vector>
#include<set>
#include <map>

using namespace std;

class scHiCs {
public:
    scHiCs(vector<string> list_of_files, string reference_genome, int resolution, bool
    adjust_resolution = true, bool sparse = false, string chromosomes = "all", string
           format = "customized", int keep_n_strata = 10, bool store_full_map = false,
           vector<string> operations = {}, int header = 0, int customized_format = 0,
           double map_filter = 0., bool gzip = false,
           bool parallelize = false, int n_processes = 0

    ) {
        this->resolution = resolution;
        this->chromosomes,
                this->chromosome_lengths = get_chromosome_lengths(reference_genome,
                                                                  chromosomes,
                                                                  resolution);
        this->num_of_cells = list_of_files.size();
        this->sparse = sparse;
        this->keep_n_strata = keep_n_strata;
        this->contacts = xt::zeros < int > {(num_of_cells)};
//        if (keep_n_strata) {
//            for (string ch : this.chromosomes) {
//                this.strata[ch];
//            }
//        }

        cout<<"Loading HiC data..."
        if(parallelize){
            throw "Not implemented yet";
        }
        else{
            int idx =0;
            for(string file : this->files){
                for(string ch : this->chromosomes){
                    size_t index = ch.find("ch");
                    if(index!=string::npos && ch.find("chr")==string::npos){
                        ch.replace(index,2,"chr");
                    }
                    xt::xarray<double>mat;
                    vector<xt::xarry<double>> strata;

                    this->contacts[idx]+=xt::sum(mat) / 2+xt::trace(mat)/2;

                }
            }

        }
    }


private:
    vector<string> files, operations;
    string reference_genome, format;
    set<string> chromosomes;
    xt::xarray<int> contacts;
    xt::xarray<double> short_range, mitotic;
    map<string, xt::xarray <xt::xarry<double>> strata;

    int
            resolution, keep_n_strata, header, customized_format, n_processes,
            chromosome_lengths, num_of_cells;
    bool
            adjust_resolution, sparse, store_full_map, gzip, parallelize;
};

#endif //SCHICTOOLS_CONTACTMAPS_H
