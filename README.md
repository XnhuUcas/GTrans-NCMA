# Ensemble graph-transfer learning under Gaussian graphical models on heterogeneous network-linked data

This repository is the official implementation of "*Ensemble graph-transfer learning under Gaussian graphical models on heterogeneous network-linked data*".

**File Structure:**

- **`data_preprocessing.R`**  
 Processes academic publication data through three steps and builds final author-term matrices and coauthor network matrices for term graph construction.

- **`method_functions.R`**  
  Defines core functions for:
  - Data generation;
  - Precision matrix estimation using GTrans-NCMA.

- **`run_example.R`**  
  Demonstrates complete workflow:
  1. Data generation;
  2. Precision matrix estimation with GTrans-NCMA;
  3. Performance evaluation using:
     - KL divergence;
     - Frobenius norm error.


## Requirements

First of all, make sure you have installed the R language environment (it is recommended to use R version 4.1 or higher).

In R, run the following command to install the required packages:

```r
install.packages(c("glasso", "Matrix", "igraph"))
```
For data preprocessing, we suggest downloading the MADStat dataset and the file 'DataForGNC-Plot-Combined.Rda' contains the top 300 terms selected by tf-idf scores, as described in the work by Li et al. (2020) on high-dimensional Gaussian graphical models for network-linked data.  The MADStat dataset is publicly available at https://github.com/ZhengTracyKe/MADStat. The file 'DataForGNC-Plot-Combined.Rda' is publicly available in the GNC repository at https://github.com/tianxili/GNC/blob/master/GNC-lasso.R.
## Training

To train the model in the paper, run the entire **`run_example.R`** script. 

## Evaluation

To evaluate my model, run:

```eval
KL_divergence <- -log(det(Omega_hat)) + sum(diag(Omega_hat %*% Sigma_true)) - (-log(det(Omega_true)) + p)
Frob_norm_error <- sum((Omega_hat - Omega_true)^2) / p  
```

## Results

The final datasets obtained after preprocessing by running **`data_preprocessing.R`** include:

| target journals  | source journals |
|------------------------------- | ----------------------------- |
|     1090 authors and 300 terms         |     8616 authors and 300 terms     |
