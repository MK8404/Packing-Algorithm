#include <algorithm>
#include <chrono>
#include <climits>
#include <iostream>
#include <queue>
#include <random>
#include <utility>
#include <vector>

#include "helper.h"

using std::cout;
using std::endl;
using std::get;
using std::make_pair;
using std::max;
using std::min;
using std::mt19937; // Standard mersenne_twister_engine seeded with rd()
using std::pair;
using std::priority_queue;
using std::random_device;
using std::sort;
using std::uniform_int_distribution;
using std::vector;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

sparse_matrix initial_solution(vector<vector<int>> A);

void shake_1(int k, sparse_matrix &spm);
void shake_2(int k, sparse_matrix &spm);
void shake(int k, int k_min, int k1_max, int k_step, sparse_matrix &spm);
void local_search(sparse_matrix &spm);
bool Move(const sparse_matrix &spm, const sparse_matrix &spm1, int alpha);

pair<vector<int>, int> VNS_Band(vector<vector<int>> A, int k_min, int k_max,
                                int k1_max, int k_step, int t_max, int alpha);
