# Packing Algorithm R Scripts

This folder contains R scripts for generating test matrices and computing the optimal permutation using the packing algorithm.

## Files

1. **`matrix_generation.R`** – Contains functions for generating three types of matrices:

   - Band matrices
   - Block tridiagonal matrices
   - Dyadic matrices
2. **`packing_alg.R`** – Implements functions necessary for computing the optimal permutation using the packing algorithm.

## Installation and Dependencies

To run the code successfully, the following R packages must be installed:

```r
install.packages(c("data.table", "Rcpp"))
```

The R-package **`data.table`** is used to speed up the computation of the alignment process for large dimensions.

The C++ interface R-package **`Rcpp`** has been used to enable C++ code for two functions:

1. **`dist_mat_rcpp`** computes the distances matrix for a given input row neighborhood (D_i) and is used in **`find_per_final`**.
2. **`neighborhood_s_cpp`** computes higher order neighborhoods and is used in **`power`**.

Utilizing C++ improves computational efficiency in these two functions.

Additionally, in function **`find_per_final`**, the multidimensional scaling function **`cmdscale`** is used. It is a part of the standard **`stats`** R-package.
