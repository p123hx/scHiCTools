//
// Created by Bean Juice on 16/12/2020.
//

#include "load_hic_file.h"
#include <cmath>
#include <set>
#include <map>
#include <string>
#include <fstream>
#include <iostream>
#include <list>
#include "processing_utils.h"
using namespace std;

static string my_path = "";

int _res_change(int length, int res) {
    return ceil(length / res);
}

struct Flg {
    vector<int> p1s, p2s;
    vector<double> vs;
};

Flg file_line_generator(string file,
                        vector<int> format =
                        vector<int>(), string
                        chrom = "", int header = 0, int resolution = 1,
                        bool resolution_adjust = true, double
                        mapping_filter = 0., bool gzip = false) {
    if (gzip) throw "gzip opening not implemented yet";
    ifstream fin;
    fin.open(file);
    Flg ans;
    if (fin.is_open()) {
        if (header) {
            for (int tmp = 0; tmp < header; tmp++) {
                string dum;
                getline(fin, dum);
            }
        }
        string line;

        while (getline(fin, line)) {
            istringstream iss(line);
            vector<std::string> lst((istream_iterator<string>(iss)),
                                    istream_iterator<string>());
            if (lst.size() != 4 && lst.size() != 5 && lst.size() != 6) {
                throw "Wrong custom format!";
            }
            if (format[0] != 0 && format[2] != 0) {
                string c1, c2 = lst[format[0] - 1],
                        lst[format[2] - 1];
                if ((c1 != chrom && "chr" + c1 != chrom) || (c2 != chrom && "chr" + c2 !=
                                                                            chrom))
                    continue;

            }
            if (format.size() == 6) {
                double q1 = stod(lst[format[4] - 1]), q2 = stod(lst[format[4] - 1]);
            }
            int p1 = stoi(lst[format[4] - 1]), p2 = stoi(lst[format[4] - 1]);
            if (resolution_adjust) {
                p1 = p1 / resolution;
                p2 = p2 / resolution;
            }
            double v;
            if (format.size() == 4 || format.size() == 6) v = 1.0;
            else if (format.size() == 5) v = stod(lst[format[4] - 1]);
            else {
                throw "Wrong custom format!";
            }
            ans.p1s.push_back(p1);
            ans.p2s.push_back(p2);
            ans.vs.push_back(v);
        }
        fin.close();

    }
    return ans;

}

pair<set<string>, map<string, int>>
get_chromosome_lengths(const string &ref_str, const string &chromosomes,
                       set<string> chrom_set, int res = 1) {
//path is in the same folder

    map<string, int> length;
    ifstream fin;
    string filename = "reference_genome/" + ref_str;
    try {
        fin.open(filename);
    } catch (const ifstream::failure &e) {
        try {
            fin.open("reference_genome\\" + ref_str);
        } catch (const ifstream::failure &e) {
            cout << "Exception opening/reading file";
        }
    }
    string name, count;
    set<string> chroms;
    while (fin >> name >> count) {
        int res_change = _res_change(stoi(count), res);
        length.insert(pair<string, int>(name, res_change));
        chroms.insert(name);
    }
    if (chromosomes == "except X") length.erase("chrX");
    else if (chromosomes == "except Y") length.erase("chrY");
    else if (chromosomes == "except XY") {
        length.erase("chrX");
        length.erase("chrY");
    }
    map<string, int> length2;

    if (!chrom_set.empty()) {
        for (string elm : chrom_set) {
            length2[elm] = length[elm];
        }
        return pair<set<string>, map<string, int >>(chrom_set, length2);
    }
    return pair<set<string>, map<string, int >>(chroms, length2);
//should we set up threshold for bad reference genome type?
}


//ref_map: char int or string int
pair<set<string>, map<string, int>>
get_chromosome_lengths(map<string, int> ref_map, string chromosomes,
                       set<string> chrom_set, int res = 1) {
    map<string, int> length;
    set<string> chroms;
    for (auto it = ref_map.begin(); it != ref_map.end(); it++) {
        length.insert(pair<string, int>(it->first, _res_change(it->second, res)));
        chroms.insert(it->first);
    }
    if (chromosomes == "except X") length.erase("chrX");
    else if (chromosomes == "except Y") length.erase("chrY");
    else if (chromosomes == "except XY") {
        length.erase("chrX");
        length.erase("chrY");
    }
    map<string, int> length2;

    if (!chrom_set.empty()) {
        for (string elm : chrom_set) {
            length2[elm] = length[elm];
        }
        return pair<set<string>, map<string, int >>(chrom_set, length2);
    }
    return pair<set<string>, map<string, int >>(chroms, length2);
}

pair<xt::xarray<double>, vector<xt::xarray<double>>> load_HiC(string file, map<string,
        int> genome_length,
                                                              string
                                                              format,
                                                              vector<int> custom_format,
                                                              int header,
                                                              string chromosome,
                                                              int resolution,
                                                              bool resolution_adjust,
                                                              double map_filter,
                                                              bool sparse,
                                                              bool gzip,
                                                              bool keep_n_strata,
                                                              vector<string>
                                                              operations) {
    int size = genome_length[chromosome];
    transform(format.begin(), format.end(), format.begin(), ::tolower);
    Flg gen;
    if (format == "shortest_score") {
        gen = file_line_generator(file,vector<int>{1,2,3,4},chromosome,0,resolution,
        resolution_adjust, 0,gzip);
    } else {
        throw "Not implemented yet";
    }
    xt::xarray<double> mat = xt::zeros<double>({size,size});
    int count = 0,gen_size = gen.vs.size();
    int p1=gen.p1s[0],p2=gen.p2s[0];double val = gen.vs[0];
    mat(p1,p2)+=val;
    if(p1!=p2) mat(p2,p2)+=val;
    for(count;count<gen_size;count++){
         p1 = gen.p1s[count]-gen.p1s[count-1];
        p2 = gen.p2s[count] - gen.p2s[count-1];
        val = gen.vs[count];
        if(count%100000==0) cout<<"Line: "<<count<<endl;
        mat(p1,p2)+=val;
        if(p1!=p2) mat(p2,p2)+=val;
    }

    if(!operations.empty())mat=matrix_operation(mat,operations);
    vector<xt::xarray<double>> strata = vector<xt::xarray<double>>();
    if(keep_n_strata){
        int matSize = mat.size();
        for(int i=0;i<keep_n_strata;i++){
            strata.push_back(xt::diag(xt::view(mat,xt::range(i,xt::placeholders::_),
                                               xt::range(xt::placeholders::_,
                                                       matSize-i))));
        }
    }
if(sparse) throw "Not implemented yet";
return make_pair(mat,strata);

}


