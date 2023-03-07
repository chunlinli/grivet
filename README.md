# Graphical Instrumental Variables Estimation and Testing
This repository contains an implementation of the following paper 

- Discovery and inference of a causal network with hidden confounding. Submitted.

The method is named **Gr**aphical **I**nstrumental **V**ariables **E**stimation and **T**esting(GrIVET).

## Contents

-  `./grivet/`: contains an R package for algorithm implementations.

- `./real_data/`: contains ADNI data and its analysis with GrIVET.

- `./simulations/*.R`: contains R scripts for simulations and summaries.

- `./simulations/primary_results`: stores simulation results.

- `./simulations/summary_results`: stores summaries of simulation results.

## Package Installations

For R, the testing version is 4.1.0. To install the package, run the following Bash script to build and install the R package.

```bash
R CMD build grivet
R CMD install grivet_1.0.tar.gz
```

Note: replace "grivet_1.0.tar.gz" by the file built if applicable.

The following packages are used.

```r
pkg <- c(
    "stats", "cvTools", " mnormt", # for R package
    "ggcorrplot", # for real data analysis
    "doParallel", "mvtnorm", "MASS", "lrpsadmm", "bnlearn","pcalg","clusterGeneration",
    "ggplot2" # for simulations
)
install.packages(pkg)
```

NOTE: some packages have dependencies unavailable from CRAN. The user may need to install them manually.

## Reproducing Simulations.

To reproduce the simulation results in the paper, change the working directory using the following Bash script.

```bash
cd simulations
```

### System information 

The code is tested on a server with specs:
```
System Version:             Ubuntu 18.04.6 LTS 4.15.0-176-generic x86_64
Model name:                 Intel(R) Xeon(R) Gold 5218 CPU @ 2.30GHz
Total Number of Cores:      64
Memory:                     528 GB
```
No GPU is required.

### Table 1

```bash
Rscript comparison1_setting1_simu.R
Rscript summary_part1_setting1.R
```
The simulation results will be stored in `./simulations/primary_results/part1/setting1` and the summary will be printed and stored in `./simulations/summary_results/part1/setting1`.








