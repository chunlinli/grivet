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

## Simulation Reproductions

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
To reproduce Table 1, run the following bash script.
```bash
Rscript comparison1_setting1_simu.R
Rscript summary_part1_setting1.R
```
The simulation results will be stored in `./simulations/primary_results/part1/setting1/` and the summary will be printed and stored in `./simulations/summary_results/part1/setting1/`.

### Table 2
To reproduce Table 2, run the following bash script.
```bash
Rscript comparison2_setting1_simu.R
Rscript summary_part2_setting1.R
```
The simulation results will be stored in `./simulations/primary_results/part2/setting1/` and the summary will be printed and stored in `./simulations/summary_results/part2/setting1/`.

## References

[1] Frot, B., Nandy, P., & Maathuis, M. H.  (2019).
[Robust causal structure learning with some hidden variables](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssb.12315), JRSSB. 
Open-sourced softwares: LRpS+GES is implemented by [lrpsadmm](https://github.com/benjaminfrot/lrpsadmm) and [pcalg](https://github.com/cran/pcalg).

[2] Colombo, D., Maathuis, M. H., Kalisch, M., & Richardson, T. S. (2012).
[Learning high-dimensional directed acyclic graphs with latent and selection variables](https://projecteuclid.org/journals/annals-of-statistics/volume-40/issue-1/Learning-high-dimensional-directed-acyclic-graphs-with-latent-and-selection/10.1214/11-AOS940.full), AOS. 
Open-sourced software: RFCI is implemented by [pcalg](https://github.com/cran/pcalg).

[3] Kalisch, M., Mächler, M., Colombo, D., Maathuis, M. H., & Bühlmann, P. (2012).
[Causal Inference Using Graphical Models with the R Package pcalg](https://www.jstatsoft.org/article/view/v047i11), JSS. 
Open-sourced software: [pcalg](https://github.com/cran/pcalg).


In addition, part of the simulation code is adapted from 
[Frot's code](https://github.com/benjaminfrot/lrpsadmm-examples)

**I would like to thank the authors of above open-sourced software.**






