# Graphical Instrumental Variables Estimation and Testing
This repository contains an implementation of the following paper 

- Discovery and inference of a causal network with hidden confounding. Submitted.

The method is named **Gr**aphical **I**nstrumental **V**ariables **E**stimation and **T**esting(GrIVET).

## Contents

-  `./grivet/`: contains an R package for algorithm implementations.

- `./real_data/`: contains analysis of ADNI data with GrIVET.

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

The simulations are tested with a request of 2 nodes and 60 GB memory on a server with specs of each node:
```
System Version:             Linux cn0203 3.10.0-1160.83.1.el7.x86_64
Model name:                 Intel(R) Xeon(R) CPU E5-2680 v3 @ 2.50GHz
CPU(s):                     24
```
No GPU is required.

### Table 1
To reproduce Table 1, run the following bash script(Approximately 34 hours).
```bash
Rscript comparison1_setting1_simu.R
Rscript summary_part1_setting1.R
```
The simulation results will be stored in `./simulations/primary_results/part1/setting1/` and the summary will be printed and stored in `./simulations/summary_results/part1/setting1/`.

### Table 2
To reproduce Table 2, run the following bash script(Approximately 21.5 hours).
```bash
Rscript comparison2_setting1_simu.R
Rscript summary_part2_setting1.R
```
The simulation results will be stored in `./simulations/primary_results/part2/setting1/` and the summary will be printed and stored in `./simulations/summary_results/part2/setting1/`.

### Table 3
To reproduce Hub Graph part of Table 3, run the following bash script(Approximately 24 hours).
```bash
Rscript infer_setting1_hub.R
Rscript summary_infer_setting1_hub.R
```
The simulation results will be stored in `./simulations/primary_results/part3/setting1/` and the summary will be printed and stored in `./simulations/summary_results/part3/setting1/`. In addition, the statistics distribution figures will be stored in `./simulations/summary_results/part3/setting1/figures/`.

To reproduce Random Graph part of Table 3, run the following bash script(Approximately 31 hours).
```bash
Rscript infer_setting1_random.R
Rscript summary_infer_setting1_random.R
```
The simulation results will be stored in `./simulations/primary_results/part3/setting1/` and the summary will be printed and stored in `./simulations/summary_results/part3/setting1/`. In addition, the statistics distribution figures will be stored in `./simulations/summary_results/part3/setting1/figures/`.

### Figure 

To reproduce Figure .., run the following bash script(Approximately 17 hours).
```bash
Rscript comparison4_discovery.R
Rscript comparison4_estimation.R
Rscript summary_part4.R
```
The simulation results will be stored in `./simulations/primary_results/part4/discovery/` and `./simulations/primary_results/part4/estimation/` and the summary of performance of structure learning will be printed and the figures will be stored in `./simulations/summary_results/part4/`.

### Figure 

To reproduce Figure .., run the following bash script(Approximately 17.5 hours).
```bash
Rscript infer_part5.R
Rscript summary_part5.R
```
The simulation results will be stored in `./simulations/primary_results/part5/` and the boxplots will be stored in  `./simulations/summary_results/part5/`.

### Table 1(Supplementary Materials)
To reproduce Table 1(Supplementary Materials), run the following bash script(Approximately 12.5 hours).
```bash
Rscript comparison1_setting2_simu.R
Rscript summary_part1_setting2.R
```
The simulation results will be stored in `./simulations/primary_results/part1/setting2/` and the summary will be printed and stored in `./simulations/summary_results/part1/setting2/`.

### Table 2(Supplementary Materials)
To reproduce Table 2(Supplementary Materials), run the following bash script(Approximately 3.5 hours).
```bash
Rscript comparison2_setting2_simu.R
Rscript summary_part2_setting2.R
```
The simulation results will be stored in `./simulations/primary_results/part2/setting2/` and the summary will be printed and stored in `./simulations/summary_results/part2/setting2/`.

### Table 3(Supplementary Materials)
To reproduce Table 3(Supplementary Materials), run the following bash script(Approximately 3.5 hours).
```bash
Rscript comparison2_setting3_simu.R
Rscript summary_part2_setting3.R
```
The simulation results will be stored in `./simulations/primary_results/part2/setting3/` and the summary will be printed and stored in `./simulations/summary_results/part2/setting3/`.

### Table 4(Supplementary Materials)
To reproduce Hub Graph part of Table 4(Supplementary Materials), run the following bash script(Approximately 2.5 hours).
```bash
Rscript infer_setting2_hub.R
Rscript summary_infer_setting2_hub.R
```
The simulation results will be stored in `./simulations/primary_results/part3/setting2/` and the summary will be printed and stored in `./simulations/summary_results/part3/setting2/`. In addition, the statistics distribution figures will be stored in `./simulations/summary_results/part3/setting2/figures/`.

To reproduce Random Graph part of Table 3(Supplementary Materials), run the following bash script(Approximately 5.5 hours).
```bash
Rscript infer_setting2_random.R
Rscript summary_infer_setting2_random.R
```
The simulation results will be stored in `./simulations/primary_results/part3/setting2/` and the summary will be printed and stored in `./simulations/summary_results/part3/setting2/`. In addition, the statistics distribution figures will be stored in `./simulations/summary_results/part3/setting2/figures/`.


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
[Frot's code](https://github.com/benjaminfrot/lrpsadmm-examples).

**I would like to thank the authors of above open-sourced software.**






