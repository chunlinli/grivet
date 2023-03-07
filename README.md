# GrIVET: Graphical Instrumental Variables Estimation and Testing
This repository contains an implementation of the following paper 

- Causal discovery and inference with confounders. Submitted.

The method is named **Gr**aphical **I**nstrumental **V**ariables **E**stimation and **T**esting(GrIVET).

## Contents

- `intdagdiscovery.R`: implements peeling algorithm for discovery of ancestral and interventional relationships.

- `intdagcoef.R`: implements coefficient estimation given the ancestral and interventional relationships.

- `mbtlp.R`: estimates precision matrix.

- `intdaginfer.r`: conducts hypothesis testing given the ancestral and interventional relationships along with a precision matrix.

- `./real_data/`: ADNI data and its analysis with GrIVET.

- `./simulations/`: simulations and simulated data.

## Preliminaries

To use the code, you need to compile C code `ctlpreg.c`. In terminal,
```
R CMD SHLIB ctlpreg.c
```

## Usage

An illustration of usage is provided in `ill_usage.Rmd`