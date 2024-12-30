# TCGAimmunosurv
# Project Overview
TCGAimmunosurv, which combines mutation-specific survival and single-cell trajectory analyses.

## Required Libraries
To ensure the package works correctly, install the following libraries (if your environment doesn't consist the following libraries):

## Installation using source:
```
R Console
source(Install.R)

```

## Data Requirements
The project requires two main datasets:

- TCGA Cancer Datasets
- Single-Cell RNA-Seq Datasets of different cancer types
-  Please find the data here: [Link](https://zenodo.org/uploads/14575108)

Note: User can download the in working directory or set the path for the data. for Bulk RNA-Seq analysis if you dont download the data it automatically start downloading from the GDC.

  ## File Descriptions and Execution

1. functions.R

Run this script after downloading data.

Execution: Run the following command.
```
R Console

source(functions.R)
```

2. eample_breast_cd8tcell.R

This file contains the eecution the functions for the analysis.

Execution: Run the commands sequential wise to reproduce the results and User can also use for the other cancer types.
