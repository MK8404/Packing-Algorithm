# **Supplementary Code for "Structural packing, factorization, and efficient inversion of sparse positive definite matrices"** by Kos, M., Podgórski, K., Wu, H.

This repository contains R scripts that serve as a supplement to the article **"Structural packing, factorization, and efficient inversion of sparse
positive definite matrices"**. The provided code enables reproduction of the simulation and computational results presented in the section dedicated to the **packing algorithm**.

## **Repository Contents**

The R Markdown files

1. **`Sect5robust.Rmd`** - the robustness study of Section 5 in the paper.
2. **`Sect5band.Rmd`** - Code needed to generate Figures in the subsection **"Permuted band matrices"**
3. **`Sect5tridiagonal.Rmd`** - Code needed to generate Figures in the subsection **"Permuted tridiagonal matrices"**
4. **`Sect5dyadic.Rmd`** - Code needed to generate Figures in the subsection **"Permuted dyadic matrices"**

The repository also includes three subdirectories described below.

### Packing

This folder contains R scripts for matrix generation and for computing the optimal permutation using the packing algorithm. For details about the scripts, installation, and dependencies, please refer to [Packing/README.md](Packing/README.md).

### VNS

This folder contains C/C++ replication code for the Variable Neighbourhood Search (VNS) algorithm for matrix bandwidth minimization described by Mladenovic et al. (2010). For details about the installation and use of the files in this folder, please refer to [VNS/README.md](VNS/README.md).

### Misc

This folder contains miscellaneous scripts to support R Markdown scripts.

## **Usage Instructions**

To reproduce the figures and simulations presented in Section 5 of the paper, please run the code contained in the R Markdown files.

---

For any questions or issues, please refer to the article or open an issue in this repository.
