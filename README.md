### **README.md**  

# **Supplementary Code for "blabla"**  

This repository contains R scripts that serve as a supplement to the article **"blabla"**. The provided code enables the reproduction of the results presented in the section dedicated to the **packing algorithm**.  

## **Repository Contents**  

The repository includes the following three R scripts:  

1. **`matrix_generation.R`** – Contains functions for generating three types of matrices:  
   - Band matrices  
   - Block tridiagonal matrices  
   - Dyadic matrices  

2. **`packing_alg.R`** – Implements functions necessary for computing the optimal permutation using the packing algorithm.  

3. **`simulations_and_plots.R`** – Provides scripts for reproducing the simulation results and generating the figures from the **"blabla"** section of the article.  

## **Installation and Dependencies**  

To run the code successfully, the following R packages must be installed:  

```r
install.packages(c("data.table", "Rcpp"))
```
The R-package **`data.table`** is used to speed up the computation of the alignment process for large dimensions.

The C++ interface R-package **`Rcpp`** has been used to enable C++ code for two functions:

**`dist_mat_rcpp`** that computes the distances matrix for a given input row neighborhood and is used in **`find_per_final`**

**`neighborhood_s_cpp`**  that computes higher order neighborhoods and is used in  **`power`**. 

Utilizing C++ improves computational efficiency in these two functions. 

Additionally, in **`find_per_final`**, the multidimensional scaling function **`cmdscale`** is used. It is a part of the standard **`stats`** R-package. 

## **Usage Instructions**  

To reproduce the results, first, ensure that all necessary functions from `matrix_generation.R` and `packing_alg.R` are loaded into the R environment. Then, execute `simulations_and_plots.R` to generate the simulation results and figures.  

---

For any questions or issues, please refer to the article or open an issue in this repository.
