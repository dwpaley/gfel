# gfel

An attempt to determine unit cells from XFEL data using PXRD indexing tools in GSAS-II.

## Overview

* GSAS-II indexing uses a modification of the SVD-Index method here: https://doi.org/10.1107/S0021889802019878. Random unit cells are generated and refined against a peak list using singular value decomposition. Figures of merit for the candidate solutions are primarily based on (1) the de Wolff M20 statistic, measuring agreement of 2th positions; and (2) the number of unindexed peaks. This search is typically run on a list of ~20 peaks.
* We have a huge number of measured d-spacings, but they are not very accurate, primarily due to the XFEL bandwidth. Our idea is to randomly sample the measured d-spacings and index the trial sets. Since the measured d-spacings center on their true values, we hope the correct unit cell will rise to the top if many random trial sets are indexed. The script `mc_cellsearch.py` takes a pickled list of reflections with noise and indexes random subsets.
* In a preliminary trial with 0.5% Gaussian noise in the d-spacings, the correct unit cell was usually in the list of candidate cells, but was not especially prominent. My next step is to take this short list of candidate cells and use them as seeds for further unit cell refinement and ranking.

## Installation

This repo contains a lightly modified version of GSASII. The main product here is the scripts (currently all in test_scripts) that import and use GSASII functions, but there are a few optimizations of the GSASII code for our purposes.

The GSASII package and its dependencies are available via conda. The idea is to put this repo in modules, use libtbx.conda to install the dependencies, and add the repo to your python import path using a .pth link in site-packages.
```
$ cd modules
$ git clone https://github.com/dwpaley/gfel
$ libtbx.conda install gsas2pkg -c briantoby --only-deps
$ echo $PWD/gfel/GSASII > ../conda_base/lib/python3.6/site-packages/GSASII.pth # modify to match your site-packages dir
```
