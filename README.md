# Soil-pH-denit
This repository contains the experimental data of soil microbiome denitrification under pH perturbations and a mathematical model of the metabolite dynamics.

doi: [https://doi.org/10.1101/2024.03.15.584851](https://doi.org/10.1101/2024.03.15.584851)

## Data
### SourceData
This folder includes the source datasets of all metabolite dynamics and sequencing data.
### ProcessedData
This folder includes the processed data that is saved in the form of MAT-file (.mat). 

The data can be read by MATLAB (r2023b or later version). 

The data is used for analysis by the MATLAB code in the Method folder.

## Method
### FunctionAnalysis
This folder includes all the codes and results for the least-square-fitting method on the metabolite dataset.

Use of the codes: fit the consumer-resource model to the nitrate time-series data. The fitted parameters are $x_0$, $C_0$, and $A_0$.

The codes can be run by MATLAB (r2023b or later version).

### AbundanceAnalysis
This folder includes all the codes and results for the coarse-graining method on the sequencing dataset.

Use of the codes: group the sequencing data at the phylum level, compute the growth folds and survival folds, then do NMF on the growth/survival matrix.

The codes can be run by MATLAB (r2023b or later version).

