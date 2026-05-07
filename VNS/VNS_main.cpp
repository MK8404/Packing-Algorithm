#include <stdio.h>
#include <stdlib.h>

#include <iostream>

#include "VNS.h"
#include "helper.h"

using std::cout;
using std::endl;
using std::string;

#include <string>

using std::stoi;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr,
            "Usage: %s [matrix-market-filename] [-k_min 5] [-k_max 100] "
            "[-k1_max 50] "
            "[-k_step 5] [-t_max 10] [-alpha 50]\n",
            argv[0]);
    exit(1);
  }

  // Default parameters matching the original Mladenovic implementation
  int k_min = 5;
  int k_max = 100;
  int k1_max = 50;
  int k_step = 5;
  int t_max = 10;
  int alpha = 50;

  for (int i = 2; i < argc; ++i) {
    string arg = argv[i];
    if (arg == "-k_min" || arg == "--k_min") {
      if (i + 1 < argc)
        k_min = stoi(argv[++i]);
    } else if (arg == "-k_max" || arg == "--k_max") {
      if (i + 1 < argc)
        k_max = stoi(argv[++i]);
    } else if (arg == "-k1_max" || arg == "--k1_max") {
      if (i + 1 < argc)
        k1_max = stoi(argv[++i]);
    } else if (arg == "-k_step" || arg == "--k_step") {
      if (i + 1 < argc)
        k_step = stoi(argv[++i]);
    } else if (arg == "-t_max" || arg == "--t_max") {
      if (i + 1 < argc)
        t_max = stoi(argv[++i]);
    } else if (arg == "-alpha" || arg == "--alpha") {
      if (i + 1 < argc)
        alpha = stoi(argv[++i]);
    }
  }

  pair<vector<vector<int>>, vector<vector<double>>> read_mm_result =
      read_mm(argc, argv);
  vector<vector<int>> A = get<0>(read_mm_result);
  vector<vector<double>> VAL = get<1>(read_mm_result);

  cout << "Running VNS with parameters: K_MIN=" << k_min << " K_MAX=" << k_max
       << " K1_MAX=" << k1_max << " K_STEP=" << k_step << " T_MAX=" << t_max
       << " ALPHA=" << alpha << endl;

  auto VNS_pair = VNS_Band(A, k_min, k_max, k1_max, k_step, t_max, alpha);
  vector<int> f_star = get<0>(VNS_pair);
  write_mm(A, VAL, f_star, argv);

  return 0;
}
