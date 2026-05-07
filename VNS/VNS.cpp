#include "VNS.h"

random_device rd;
mt19937 gen(rd());

sparse_matrix initial_solution(vector<vector<int>> A) {
  // f should be of length n.
  int n = A.size();
  vector<int> f(n, -1);
  vector<bool> mark(n, false);

  // the random number generator
  uniform_int_distribution<> runif(0, n - 1);
  int k = runif(gen);
  mark[k] = true;

  // initialization
  vector<int> q;
  q.emplace_back(k);
  int ql = 1;
  f[k] = 0;
  int l = 0;
  int ptr = 0;

  // depth first search
  while (l < n) {
    // cout << l << endl;
    int r = 0;
    vector<int> s;
    bool is_new = false;
    for (int i1 = 0; i1 < ql; i1++) {
      int i = q[i1];
      for (int j1 = 0; j1 < A[i].size(); j1++) {
        int j = A[i][j1];
        if (not mark[j]) {
          r++;
          s.emplace_back(j);
          mark[j] = true;
          is_new = true;
        }
      }
    }

    uniform_int_distribution<> runif_temp(0, ql - 1);
    int j_star = runif_temp(gen);

    for (int j1 = 0; j1 < ql; j1++) {
      int j = q[j_star];
      f[j] = l;
      l = l + 1;
      j_star = j_star + 1;
      if (j_star > ql - 1) {
        j_star = 0;
      }
    }

    if (!is_new) {
      for (; ptr < n; ptr++) {
        if (!mark[ptr]) {
          s = vector<int>({ptr});
          r = 1;
          mark[ptr] = true;
          break;
        }
      }
    }

    q = s;
    ql = r;
  }

  return sparse_matrix(A, f);
}

template <class T, class S, class C> S &Container(priority_queue<T, S, C> &q) {
  struct HackedQueue : private priority_queue<T, S, C> {
    static S &Container(priority_queue<T, S, C> &q) {
      return q.*&HackedQueue::c;
    }
  };
  return HackedQueue::Container(q);
}

void shake_1(int k, sparse_matrix &spm) {
  // find the corresponding set K
  priority_queue<pair<int, int>, vector<pair<int, int>>, decltype(cmp_pair) *>
      K(cmp_pair);
  int n = spm.n;
  for (int i = 0; i < n; i++) {
    pair<int, int> t_i = make_pair(i, get<0>(spm.B_f[i]));
    if (K.size() < k) {
      K.push(t_i);
    } else {
      if (!cmp_pair(K.top(), t_i)) {
        K.pop();
        K.push(t_i);
      }
    }
  }

  vector<pair<int, int>> vec_K(Container(K));

  // K is determined. Now we begin to swap.
  for (int i = 0; i < k; i++) {
    uniform_int_distribution<> runif_1(0, vec_K.size() - 1);
    int k = runif_1(gen);
    pair<int, int> tup_u = vec_K[k];
    vec_K.erase(vec_K.begin() + k);
    int u = get<0>(tup_u);

    // find v (without randomness for multiple solutions)
    int v = get<1>(spm.B_f[u]);

    // find w
    int w = -1;
    int min_max_wv = n;
    for (int w1 = 0; w1 < n; w1++) {
      if (spm.f_max[u] >= spm.f[w1] && spm.f_min[u] <= spm.f[w1]) {
        int max_wv = max(spm.f_max[w1] - spm.f[v], spm.f[v] - spm.f_min[w1]);
        if (max_wv < min_max_wv) {
          min_max_wv = max_wv;
          w = w1;
        }
      }
    }

    // swap the f value
    if (w != -1) {
      spm.f_swap(v, w);
    }
  }
}

void shake_2(int k, sparse_matrix &spm) {
  uniform_int_distribution<> runif_1(
      0, spm.n - 3); // so that we will always find a proper m;
  uniform_int_distribution<> runif_2(2, min(20, spm.bandwidth / 2));
  for (int i = 0; i < k; i++) {
    int b = runif_1(gen);
    int e = b + runif_2(gen);
    if (e > spm.n - 1)
      e = spm.n - 1;
    uniform_int_distribution<> runif_3(1, e - b - 1);
    int m = runif_3(gen);
    vector<int> vs(e - b + 1);
    vector<int> vals(e - b + 1);
    // #pragma omp parallel for num_threads(4)
    for (int j = b; j < e + 1; ++j) {
      vs[j - b] = spm.f_inverse[j];
      if (j >= b + m) {
        vals[j - b] = j - m;
      } else {
        vals[j - b] = j + e - b - m + 1;
      }
    }
    spm.f_assign(vs, vals);
  }
}

void shake(int k, int k_min, int k1_max, int k_step, sparse_matrix &spm) {
  if (k <= k1_max) {
    shake_1(k, spm);
  } else {
    int k1 = int((k - k_min) / k_step);
    shake_2(k1, spm);
  }
}

// use the improved hill climbing algorithm
void local_search(sparse_matrix &spm) {
  bool flag = true;
  vector<pair<int, int>> N_1;
  N_1.reserve(spm.n);
  while (flag) {
    flag = false;
    for (int v = 0; v < spm.n; ++v) {
      if (spm.critical_vertices.count(v)) {
        int f_mid = (spm.f_max[v] + spm.f_min[v]) / 2;
        int v_mid_v = abs(f_mid - spm.f[v]);
        N_1.clear();
        for (int u = 0; u < spm.n; u++) {
          int u_mid_v = abs(f_mid - spm.f[u]);
          if (u_mid_v < v_mid_v)
            N_1.push_back(make_pair(u, u_mid_v));
        }
        sort(N_1.begin(), N_1.end(),
             [](const pair<int, int> &a, const pair<int, int> &b) {
               return a.second < b.second;
             });

        for (auto const &item : N_1) {
          int u = item.first;
          int c_u = spm.critical_vertices.count(u);
          int c_v = 1;
          int c1_u = 0;
          int c1_v = 0;
          swap(spm.f[u], spm.f[v]);
          for (int w : spm.A[u]) {
            if (abs(spm.f[w] - spm.f[u]) >= spm.bandwidth) {
              c1_u = 2;
              break;
            };
            if (abs(spm.f[w] - spm.f[u]) == spm.bandwidth)
              c1_u = 1;
          }
          for (int w : spm.A[v]) {
            if (abs(spm.f[w] - spm.f[v]) > spm.bandwidth) {
              c1_v = 2;
              break;
            };
            if (abs(spm.f[w] - spm.f[v]) == spm.bandwidth)
              c1_v = 1;
          }
          swap(spm.f[u], spm.f[v]);
          if (c1_u <= c_u && c1_v <= c_v && c1_u + c1_v < c_u + c_v) {
            spm.f_swap(v, u);
            flag = true;
            break;
          }
        }
      }
    }
  }
  // auto t2 = high_resolution_clock::now();
  // auto ms_int = duration_cast<microseconds>(t2 - t1);
  // cout << "Local Search takes " << double(ms_int.count()) / 1000 << " ms"
  //      << endl;
}

bool Move(const sparse_matrix &spm, const sparse_matrix &spm1, int alpha) {
  // condition 1
  if (spm.bandwidth < spm1.bandwidth) {
    return false;
  } else if (spm.bandwidth > spm1.bandwidth) {
    // cout << "new bw = " << spm1.bandwidth << endl;
    return true;
  } else {
    // condition 2;
    if (spm1.critical_vertices.size() < spm.critical_vertices.size()) {
      // cout << "new cr = " << spm1.critical_vertices.size() << endl;
      return true;
    }
    // condition 3;
    int rho = -1;
    for (int i = 0; i < spm.n; i++) {
      if (spm.f[i] != spm1.f[i]) {
        rho++;
      }
    }
    if (rho > alpha) {
      // cout << "rho = " << rho << endl;
      return true;
    } else {
      return false;
    }
  }
}

pair<vector<int>, int> VNS_Band(vector<vector<int>> A, int k_min, int k_max,
                                int k1_max, int k_step, int t_max, int alpha) {
  high_resolution_clock::time_point t1, t2;
  milliseconds ms_int;
  int B_star = INT_MAX;
  int t = 0;
  int i_max = (k_max - k_min) / k_step;
  int n = A.size();
  vector<int> f_star(n);
  for (int j = 0; j < n; ++j)
    f_star[j] = j;
  cout << "Initial bandwidth = " << bandwidth(A, f_star) << endl;
  while (t < t_max) {
    t1 = high_resolution_clock::now();
    cout << "Epoch " << t << "... " << endl;
    vector<int> f(n, -1);
    sparse_matrix spm = initial_solution(A);
    local_search(spm);
    int i = 0;
    int k = k_min;
    while (i < i_max) {
      sparse_matrix spm1(spm);
      shake(k, k_min, k1_max, k_step, spm1);
      local_search(spm1);
      if (Move(spm, spm1, alpha)) {
        spm = spm1;
        k = k_min;
        // i = 0; // Removing this to prevent infinite loops on plateaus
        i++;
      } else {
        k += k_step;
        i++;
      }
    }
    if (spm.bandwidth < B_star) {
      B_star = spm.bandwidth;
      f_star = spm.f;
    }
    cout << "SPM Bandwidth = " << spm.bandwidth << endl;
    cout << "Best Bandwidth = " << B_star << endl;
    t++;
    t2 = high_resolution_clock::now();
    ms_int = duration_cast<milliseconds>(t2 - t1);
    cout << "Epoch " << t - 1 << " takes " << ms_int.count() << " ms" << endl;
  }
  return make_pair(f_star, B_star);
}