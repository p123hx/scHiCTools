//
// Created by Bean Juice on 16/12/2020.
//

#include "load_hic_file.h"
#include <cmath>
#include <set>
#include <map>
#include <string>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <list>
using namespace std;
static string my_path = "";

int _res_change(int length, int res) {
    return ceil(length / res);
}

pair<set<string>, map<string, int>>
get_chromosome_lengths(const string &ref_str, const string &chromosomes, set<string> chrom_set,int res = 1) {
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
    set<string>chroms;
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

    if(!chrom_set.empty()){
        for(string elm : chrom_set){
            length2[elm] = length[elm];
        }
        return pair<set<string>, map<string, int>>(chrom_set,length2);
    }
    return  pair<set<string>, map<string, int>>(chroms,length2);
//should we set up threshold for bad reference genome type?
}



//ref_map: char int or string int
pair<set<string>, map<string, int>>
get_chromosome_lengths(map<string, int> ref_map, string chromosomes, set<string> chrom_set,int res = 1) {
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

    if(!chrom_set.empty()){
        for(string elm : chrom_set){
            length2[elm] = length[elm];
        }
        return pair<set<string>, map<string, int>>(chrom_set,length2);
    }
    return  pair<set<string>, map<string, int>>(chroms,length2);
}

