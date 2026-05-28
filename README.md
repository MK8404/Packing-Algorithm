### **README.md**  

# **Supplementary Code for "Structural packing, factorization, and efficient inversion of sparse positive definite matrices"** by Kos, M., Podgórski, K., Wu, H.

This repository contains R scripts that serve as a supplement to the article **"Structural packing, factorization, and efficient inversion of sparse
positive definite matrices"**. The provided code enables reproduction of the simulation and computational results presented in the section dedicated to the **packing algorithm**.  

## **Repository Contents**  
The R Markdown files
1. **`Sect5robust.Rmd`** - the robustness study of Section 5 in the paper.
2. **`Sect5band.Rmd`** - Code needed to generate Figures in the subsection **"Permuted band matrices"**
3. **`Sect5tridiagonal.Rmd`** - Code needed to generate Figures in the subsection **"Permuted tridiagonal matrices"**
4. **`Sect5dyadic.Rmd`** - Code needed to generate Figures in the subsection **"Permuted dyadic matrices"**
    

The repository also includes three subdirectories.

**'Packing'** contains 

   - R scripts:  

      1. **`matrix_generation.R`** – Contains functions for generating three types of matrices:  
         - Band matrices  
         - Block tridiagonal matrices  
         - Dyadic matrices  

      2. **`packing_alg.R`** – Implements functions necessary for computing the optimal permutation using the packing algorithm.  

**'VNS'** contains 

   - Cpp scripts with the implementation of the Variable Neighbourhood Search algorithm.  

**'Misc'** contains 

   - miscellaneous scripts to support R Markdown scripts


## **Installation and Dependencies**  

To run the code successfully, the following R packages must be installed:  

```r
install.packages(c("data.table", "Rcpp"))
```
The R-package **`data.table`** is used to speed up the computation of the alignment process for large dimensions.

The C++ interface R-package **`Rcpp`** has been used to enable C++ code for two functions:

**`dist_mat_rcpp`** that computes the distances matrix for a given input row neighborhood (D_i) and is used in **`find_per_final`**

**`neighborhood_s_cpp`**  that computes higher order neighborhoods and is used in  **`power`**. 

Utilizing C++ improves computational efficiency in these two functions. 

Additionally, in **`find_per_final`**, the multidimensional scaling function **`cmdscale`** is used. It is a part of the standard **`stats`** R-package. 

## **Usage Instructions**  

To reproduce the results, first, ensure that all necessary functions from `matrix_generation.R` and `packing_alg.R` are loaded into the R environment. Then, execute `simulations_and_plots.R` to generate the simulation results and figures.  

---

For any questions or issues, please refer to the article or open an issue in this repository.
