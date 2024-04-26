
# European

<!-- badges: start -->

[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://choosealicense.com/licenses/mit/)
<!-- badges: end -->

![](figures/figure0.jpg)

## Overview

## Data sources

![](figures/figure_data_sources.jpg)

This project uses the following databases:

| Creator | Repository | Program | Area | References |
|:--------|:-----------|:--------|:----:|:----------:|
|         |            |         |      |            |
|         |            |         |      |            |
|         |            |         |      |            |
|         |            |         |      |            |

## Workflow

The analysis pipeline follows these steps:

1.  Download eddy data
2.  Download acoustic data
3.  Overlap the two data sources to only keep acoustic data inside and
    around eddies
4.  Make and Export the Figures 1 and 2
5.  Compare inside acoustic values to outside acoustic values in each
    eddy
6.  Make and Export the Figures 3 and 4
7.  Compare the spatial distribution of the results with the previously
    published results

## Content

This repository is structured as follow:

- [`data/`](https://github.com/auroreRECE/eddy_micronecton/tree/main/data):
  contains a sampled of data used in the article. The folder is split in
  two folders:
  [`data/raw/`](https://github.com/auroreRECE/eddy_micronecton/tree/main/data/raw)
  to find a sample of all raw data (acoustic, eddy trajectories, chloro
  and sst) ; and
  [`data/intermediate/`](https://github.com/auroreRECE/eddy_micronecton/tree/main/data/intermediate)
  with some intermediate .Rdata files, to make easily the Figures.

- [`figures/`](https://github.com/auroreRECE/eddy_micronecton/tree/main/figures):
  contains the figures used to validate et visualize the outputs.

- [`scripts/`](https://github.com/auroreRECE/eddy_micronecton/tree/main/scripts):
  contains R scripts to run the workflow. The order to run these scripts
  is explained in each name of files and follow the Workflow
  description.

## Citation

Please use the following citation:

> Receveur A, Menkes C, Lengaigne M, Ariza A, Bertrand A, Dutheil C,
> Cravatte S, Allain A, Barbin L, Lebourges-Dhaussy A, Lehodey P, Nicol
> S. (2024) Code for “A rare oasis effect for forage fauna in oceanic
> eddies at the global scale.” URL:
> <https://github.com/auroreRECE/eddy_micronecton/tree/main/>.

## Contributing

All types of contributions are encouraged and valued.

## Acknowledgments

Aurore Receveur is part of the MAESTRO group co-funded by the Centre for
the Synthesis and Analysis of Biodiversity (CESAB) of the Foundation for
Research on Biodiversity (FRB), and by France Filière Pêche.
