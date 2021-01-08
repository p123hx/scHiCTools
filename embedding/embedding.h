//
// Created by Bean Juice on 06/01/2021.
//

#ifndef SCHICTOOLS_EMBEDDING_H
#define SCHICTOOLS_EMBEDDING_H

#include <cmath>
#include <set>
#include <map>
#include <string>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <list>

using namespace std;

pair<set<string>, map<string, int>>
get_chromosome_lengths(const string &ref_str, const string &chromosomes,
                       set<string> chrom_set, int res = 1);
pair <set<string>, map<string, int>>
get_chromosome_lengths(map<string, int> ref_map, string chromosomes,
                       set <string> chrom_set, int res = 1);
#endif //SCHICTOOLS_EMBEDDING_H
