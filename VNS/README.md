# Variable Neighbourhood Search for Bandwidth Minimization

This folder contains C/C++ replication code for the Variable Neighbourhood Search (VNS) algorithm for matrix bandwidth minimization described by Mladenovic et al. (2010).

The program reads a sparse matrix in Matrix Market format, applies VNS to find a bandwidth-reducing permutation, and writes the permuted matrix and permutation to Matrix Market files.

## Files

- `VNS_main.cpp`: command-line entry point
- `VNS.cpp`, `VNS.h`: VNS implementation
- `helper.cpp`, `helper.h`: Matrix Market I/O and utility routines
- `mmio.c`, `mmio.h`: Matrix Market reader/writer support

## Build

Requirements:

- A C++ compiler with C++11 support

### macOS

From this directory, compile with:

```sh
g++ -std=c++11 VNS_main.cpp VNS.cpp helper.cpp -x c++ mmio.c -o VNS_main
```

### Windows

Windows users can build the program with MSYS2/MinGW.

1. Install MSYS2: https://www.msys2.org/
2. Open the **MSYS2 MinGW 64-bit** shell.
3. Install the required packages:

```sh
pacman -S --needed mingw-w64-x86_64-gcc
```

4. Compile from this directory:

```sh
g++ -std=c++11 VNS_main.cpp VNS.cpp helper.cpp -x c++ mmio.c -o VNS_main.exe
```

Then run the executable from the same shell:

```sh
./VNS_main.exe path/to/matrix.mtx
```

## Usage

```sh
./VNS_main path/to/matrix.mtx
```

Optional parameters:

```sh
./VNS_main path/to/matrix.mtx -k_min 5 -k_max 100 -k1_max 50 -k_step 5 -t_max 10 -alpha 50
```

Defaults match the original implementation:

- `k_min = 5`
- `k_max = 100`
- `k1_max = 50`
- `k_step = 5`
- `t_max = 10`
- `alpha = 50`

## Output

For an input file named `matrix.mtx`, the program writes:

- `matrix_converted.mtx`: matrix after applying the computed permutation
- `matrix_perm.mtx`: permutation matrix

Output files are written to the current working directory.

## References

- Mladenovic, N., Urosevic, D., Perez-Brito, D., and Garcia-Gonzalez, C. G. (2010). Variable neighbourhood search for bandwidth reduction. *European Journal of Operational Research*, 200(1), 14-27. https://doi.org/10.1016/j.ejor.2008.12.015
- Matrix Market I/O library for ANSI C. http://math.nist.gov/MatrixMarket/
