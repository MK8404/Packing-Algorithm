#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <unordered_set>
#include <utility>
#include <vector>

#include "mmio.h"

using std::abs;
using std::copy;
using std::cout;
using std::endl;
using std::get;
using std::inserter;
using std::make_pair;
using std::pair;
using std::round;
using std::swap;
using std::unordered_set;
using std::vector;

#pragma once
class sparse_matrix {
   public:
    vector<vector<int>> A;
    vector<int> f;
    vector<int> f_inverse;
    vector<pair<int, int>> B_f;
    vector<int> f_min;
    vector<int> f_max;
    unordered_set<int> critical_vertices;
    int n;
    int bandwidth;

    void init_helper(vector<vector<int>> A, vector<int> f);
    sparse_matrix(const sparse_matrix &spm);
    sparse_matrix(vector<vector<int>> A);
    sparse_matrix(vector<vector<int>> A, vector<int> f);
    void flush_info(int v);
    void flush_bandwidth();
    void f_swap(int v, int u);
    void f_assign(vector<int> vs, vector<int> vals);
    ~sparse_matrix();
};

pair<vector<vector<int>>, vector<vector<double>>> read_mm(int argc,
                                                          char *argv[]);
void write_mm(vector<vector<int>> A, vector<vector<double>> VAL,
              vector<int> f_star, char *argv[]);
int bandwidth(vector<vector<int>> A, vector<int> &f);
bool cmp_pair(const pair<int, int> t1, const pair<int, int> t2);
