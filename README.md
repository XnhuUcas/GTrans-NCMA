# Ensemble graph-transfer learning under Gaussian graphical models on heterogeneous network-linked data

This repository is the official implementation of the paper:

**"Ensemble graph-transfer learning under Gaussian graphical models on heterogeneous network-linked data"**

This repository provides the complete R implementation of the proposed **GTrans-NCMA** method, including:

- The implementation of the proposed method;
- Numerical simulation experiments;
- Data preprocessing based on the MADStat publication dataset.

The proposed method aims to improve precision matrix estimation for heterogeneous Gaussian graphical models by incorporating network information and auxiliary studies through an ensemble transfer learning framework.

---

# File Structure

The repository is organized as follows:

```
.
├── README.md
│
├── run_example.R
│
├── method_functions.R
│
├── simulation_code/
│   ├── *.R
│   ├── method_functions.R
│   ├── CLIME-opt.R
│   ├── GNC-lasso.R
│   └── maegg.R
│
└── empirical_analysis_code/
    └── data_preprocessing.R
```

---

## `run_example.R`

This script demonstrates the complete workflow of the proposed method.

The workflow includes:

1. Data generation;
2. Precision matrix estimation using **GTrans-NCMA**;
3. Evaluation of estimation accuracy using:
   - Kullback–Leibler (KL) divergence;
   - Frobenius norm error.

---

## `method_functions.R`

This file contains the core functions of the proposed method.

It includes:

- Data generation functions;
- Precision matrix estimation using **GTrans-NCMA**;


---

## `simulation_code/`

This folder contains all R codes used to reproduce the numerical simulation experiments reported in the paper.

Each R script corresponds to a specific simulation setting, figure, or table in the manuscript.

This folder also contains implementations of several competing methods:

### `CLIME-opt.R`

Implementation of the CLIME estimator:

> Cai, T., Liu, W. & Luo, X. (2011).  
> A constrained ℓ1 minimization approach to sparse precision matrix estimation.  
> *Journal of the American Statistical Association*, 106(494), 594–607.

---

### `GNC-lasso.R`

Implementation of the GNC-Lasso method:

> Li, T., Qian, C., Levina, E. & Zhu, J. (2020).  
> High-dimensional Gaussian graphical models on network-linked data.  
> *Journal of Machine Learning Research*, 21(74), 1–45.

---

### `maegg.R`

Implementation of the MAEGG method:

> Liu, H. & Zhang, X. (2023).  
> Frequentist model averaging for undirected Gaussian graphical models.  
> *Biometrics*, 79(3), 2050–2062.

---

## `empirical_analysis_code/`

This folder contains the data preprocessing  code for the empirical analysis based on the MADStat publication dataset.

The main preprocessing file is:

```
data_preprocessing.R
```

This script implements the data preprocessing procedure for **Section 5.2: Trend Analysis of Research Topics** in the paper.

The preprocessing procedure includes:

1. Construction of author-term matrices;
2. Construction of coauthor network matrices;


Note that the preprocessing procedure for **Section 5.1: Performance Comparison** follows the same framework. The only difference is the selection of publication time periods. Therefore, Section 5.1 and Section 5.2 share the same preprocessing strategy, while different time windows are used to construct the corresponding datasets.

---

# Requirements

The implementation is based on the R programming language.

It is recommended to use:

```
R version >= 4.1
```

Additional packages may be required when running simulation experiments or empirical analyses.

---

# Data Availability

## MADStat Dataset

The empirical analysis uses the MADStat publication dataset.

The dataset is publicly available at:

https://github.com/ZhengTracyKe/MADStat

---

# Usage

## 1. Run Example

To reproduce the basic workflow of the proposed method, run:

```r
source("run_example.R")
```

This script demonstrates:

- Data generation;
- Application of GTrans-NCMA;
- Evaluation of estimation performance.

---

## 2. Reproduce Simulation Studies

All numerical simulation experiments in the paper can be reproduced using the scripts in:

```
simulation_code/
```

Each R script corresponds to a specific numerical experiment, figure, or table reported in the manuscript.

For example:

```r
source("simulation_code/xxx.R")
```

where `xxx.R` should be replaced by the corresponding simulation script.

---

## 3. Reproduce Empirical Analysis

The data preprocessing in **Section 5.2: Trend Analysis of Research Topics** can be reproduced using:

```
empirical_analysis_code/data_preprocessing.R
```

Running this script generates the processed author-term matrices and coauthor network matrices required for the empirical analysis.
