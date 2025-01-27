# supplementary-code-R-ClickClust of the paper:

This repository contains supplementary code for clustering categorical data using the Expectation-Maximization (EM) algorithm, implemented within the ClickClust package in R. The code accompanies the paper:

Maqsood, R., Ceravolo, P., Romero, C., Ventura, S. (2022). Modeling and predicting students’ engagement behaviors using mixture Markov models. *Knowledge and Information Systems*, 64, 1349–1384. [https://doi.org/10.1007/s10115-022-01674-9](https://doi.org/10.1007/s10115-022-01674-9)

Supplementary functions for clustering categorical data using EM algorithm implemented in ClickClust package of R.

Classifier code for mixture Markov models (implemented for ClickClust package, same file structure is used)

Following three files (containing relevant code only) are shared here: "main.R", "libClickClust.c", "libEM.c"

## Structure

- `README.md`: Provides an overview and instructions for using the code.
- `libClickClust.c`: C source code implementing specific functionalities for click pattern analysis.
- `libEM.c`: C source code implementing the Expectation-Maximization algorithm for parameter estimation.
- `main.R`: Main R script containing the workflow to load data, fit models, and evaluate results.

## Requirements

- **R**: Version 4.0 or higher.
- **C Compiler**: GCC or Clang for compiling the C source files.
- **R Packages**: Ensure the following packages are installed:
  - `ClickClust`: For clustering categorical sequences.
  - `Rcpp`: To facilitate integration of R and C++.
  - `RcppArmadillo`: For linear algebra operations.

Install the required R packages using:

```R
install.packages(c("ClickClust", "Rcpp", "RcppArmadillo"))
```

## Installation

1. **Clone the repository**:

   ```git clone https://github.com/r-maqsood/Mixture-Markov-Models-R.git
   cd Mixture-Markov-Models-R
   ```

2. **Compile the C source files**:

   Use the following commands to compile the C source files:

   ```
   R CMD SHLIB libClickClust.c
   R CMD SHLIB libEM.c
   ```


## Usage

1. **Load the compiled shared objects in R**:

   ```R
   dyn.load("libClickClust.so")  # or "libClickClust.dll" on Windows
   dyn.load("libEM.so")          # or "libEM.dll" on Windows
   ```

2. **Run the main script**:

   Open and execute `main.R` in R. This script demonstrates how to use the implemented functions for clustering categorical data using the EM algorithm.

## Reference

Maqsood, R., Ceravolo, P., Romero, C., Ventura, S. (2022). Modeling and predicting students’ engagement behaviors using mixture Markov models. *Knowledge and Information Systems*, 64, 1349–1384. [https://doi.org/10.1007/s10115-022-01674-9](https://doi.org/10.1007/s10115-022-01674-9)

```
@article{maqsood2022modeling,
  title={Modeling and predicting students’ engagement behaviors using mixture Markov models},
  author={Maqsood, Rabia and Ceravolo, Paolo and Romero, Crist{\'o}bal and Ventura, Sebasti{\'a}n},
  journal={Knowledge and Information Systems},
  volume={64},
  number={5},
  pages={1349--1384},
  year={2022},
  publisher={Springer}
}
```
