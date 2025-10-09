library(glasso)

source("method_functions.R")

K <- 5       # Number of source domains
Info.size <- 3 # Number of informative source models
p <- 20    # Dimension
n.vec <- c(100, rep(225, K))  # Target data size (100), the data size for each source (225)
Theta0.type <- 'Random' # Structure of the target precision matrix

set.seed(123)
#######################  Data generation #######################

sim.data <- DataGenerate(K = K, Info.size = Info.size, h = p * 0.4, n.vec = n.vec, p = p, type = Theta0.type)
X <- sim.data$X
Omega_true <- sim.data$Theta0 # True precision matrix of the target dataset
Sigma_true <- solve(Omega_true)

# Data split by source
X.split <- list()
start_idx <- 1
for (k in 1:(K + 1)) {
  end_idx <- start_idx + n.vec[k] - 1
  X.split[[k]] <- X[start_idx:end_idx, ]
  start_idx <- end_idx + 1
}

# Generation of mean matrices for each source
Data.list <- list()
for (k in 1:length(X.split)) {
  X.class <- X.split[[k]]
  n <- nrow(X.class)
  m <- sqrt(n)
  GridData <- MeanmatrixGen(m = m)
  Adj <- GridData$A
  D <- diag(colSums(Adj))
  L <- D - Adj
  eig <- eigen(L, symmetric = TRUE)
  m1 <- eig$vectors[, n]
  m2 <- eig$vectors[, n - 1]
  mm <- sqrt(0.95) * m1 + sqrt(0.05) * m2
  M <- mm * sqrt(n) * 1.5
  M <- matrix(rep(M, p), ncol = p)
  Data.list[[k]] <- list(X.class = X.class, GridData = GridData, M = M)
}

# Data centering
new.X.list <- list()
for (k in 1:length(X.split)) {
  X.class <- X.split[[k]]
  M <- Data.list[[k]]$M
  new.X <- X.class + M
  new.X.std <- scale(new.X, center = TRUE, scale = FALSE)
  new.X.list[[k]] <- new.X.std
}

X.all <- do.call(rbind, new.X.list)

####################### GTrans-NCMA #######################

# Potential auxiliary graph estimation
Omega.list <- list()
for (k in 1:length(new.X.list)) {
  X.std <- new.X.list[[k]]
  n <- nrow(X.std)
  Adj <- Data.list[[k]]$GridData$A
  d <- mean(colSums(Adj))
  alpha.cv <- sqrt(n)
  const <- 0.5
  best.lambda <- const * 2 * sqrt(log(p) / n)
  Omega_hat2 <- GTrans_NCMA_single(X = X.std, A = Adj, lambda = best.lambda, alpha = alpha.cv, tol = 1e-4, verbose = FALSE)
  Omega.list[[k]] <- Omega_hat2
}

# Graph-transfer by weighted combination
GTrans_NCMA <- GTrans_NCMA(X.all[1:n.vec[1], ], Data.list[[1]]$GridData$A, Omega.list, CV_num_group = 5, V = 5)
Omega_hat <- GTrans_NCMA$Omega_hat # Estimation of target precision matrix by GTrans-NCMA
Gamma_hat <- GTrans_NCMA$gamma_hat # Weights for (target, 3 informative sources, 2 non-informative sources)

# KL divergence and normalized Frobenius norm error calculation
KL_divergence <- -log(det(Omega_hat)) + sum(diag(Omega_hat %*% Sigma_true)) - (-log(det(Omega_true)) + p)# KL divergence
Frob_norm_error <- sum((Omega_hat - Omega_true)^2) / p  # Normalized Frobenius norm error

KL_divergence
Frob_norm_error