#' @title main function for maegg
#' @description Please input data and candidate models, and it'll output estimation of precision matrix
#' @param X a matrix of observations of variables.
#' @param CV_num_group number of fold used for cross-validation
#' @param glassopath_X solution path given by GLasso or similar list of graphs
#' @param num_potential_model number of potential candidate models
#' @param sample_covariance sample covariance matrix
#' @return a list including: model_set, weight, Omega_hat
#' @import pracma
#' @import huge
#' @export
maegg <- function(X, CV_num_group=5, glassopath_X=NULL, num_potential_model=NULL, sample_covariance=NULL) {

  p <- ncol(X)
  n <- nrow(X)

  num_distinct_model_in_glassopath <- 1             # the first is zero-model.
  model_space_duplicated <- array(0, dim=c(p*(p-1)/2, num_potential_model))
  for (idx_current_model in c(2:num_potential_model)) {
    current_theta_from_glasso <- glassopath_X$wi[,,idx_current_model]
    current_model <- current_theta_from_glasso[upper.tri(current_theta_from_glasso, diag = FALSE)]
    if (! identical(model_space_duplicated[, num_distinct_model_in_glassopath], current_model)) {
      model_space_duplicated[, num_distinct_model_in_glassopath+1] <- current_model
      num_distinct_model_in_glassopath <- num_distinct_model_in_glassopath + 1
    }
  }
  model_space <- model_space_duplicated[, 1:num_distinct_model_in_glassopath]
  model_space <- cbind(model_space, rep(1,p*(p-1)/2))      # add full model into the space
  num_model <- dim(model_space)[2]

  # Constrained MLE
  estimate_Omega <- function(graphical_model, covariance_train){
    # Omega is estimated by ELS Algorithm 17.1 in page 634
    p <- dim(covariance_train)[1]               # this p is the dim being chosen
    Omega_hat <- matrix(rep(0, p*p), p)
    # beta_hat is the core parameter to obtain Omega_hat
    beta_hat <- matrix(rep(0, p*p),p) # each col represent a node
    # adjacency_matrix is a 0-1 value matrix with size of p*p
    adjacency_matrix <- matrix(rep(0, p*p), p)
    adjacency_matrix[upper.tri(adjacency_matrix)] <- graphical_model
    adjacency_matrix <- adjacency_matrix + t(adjacency_matrix) + diag(p)
    # set loop
    convergence_status <- rep(FALSE, p)
    current_idx <- 0
    W_current <- W_next <- covariance_train
    iter_num <- 0
    while (!identical(convergence_status, rep(TRUE, p))) {
      current_idx <- current_idx + 1 - (current_idx %/% p)*p
      adjacency_nodes_with_current <- setdiff(which(adjacency_matrix[current_idx,] != 0), current_idx)
      # step (a): Partition W
      W11 <- W_next[-current_idx, -current_idx]
      W11_star <- W_next[adjacency_nodes_with_current, adjacency_nodes_with_current]
      s12_star <- covariance_train[current_idx, adjacency_nodes_with_current]
      if (length(adjacency_nodes_with_current) > 0){
        # it is possible that current node does not connect to any other node,
        # then W11_star is not reversible.
        # step (b)
        beta_star <- solve(as.matrix(W11_star)) %*% s12_star
        beta_hat[adjacency_nodes_with_current, current_idx] <- beta_star
      }
      # step (C), W12 is feasible, because beta_hat is zero vector.
      W12 <- W11 %*% beta_hat[-current_idx, current_idx]
      # Update W_next
      W_next[-current_idx, current_idx] <- W12
      iter_num <- iter_num + 1
      # Termination condition, converge for every idx
      if (norm(W_current-W_next, "F") < 1e-5) {
        # This condition is easy to achieve.
        convergence_status[current_idx] = TRUE
        W_hat <- W_next
      }
      W_current <- W_next
    }
    # Calculate Omega_hat, Omega12 and Omega22
    for (current_idx in c(1:p)) {
      s22 <- covariance_train[current_idx, current_idx]
      W12 <- W_hat[-current_idx, current_idx]
      Omega22 <- 1 / (s22 - W12 %*% beta_hat[-current_idx, current_idx])
      Omega12 <- - beta_hat[-current_idx,current_idx] * as.numeric(Omega22)
      # Omega12 and Omega21 and Omega22 could be updated.
      Omega_hat[-current_idx, current_idx] <- Omega12
      Omega_hat[current_idx, -current_idx] <- Omega12
      Omega_hat[current_idx, current_idx] <- Omega22
    }
    # # debug
    # print(norm(Omega_hat-solve(covariance_train), "F"))
    # # end debug
    gc()    # clean memory

    return(Omega_hat)
  }

  # # CV Procedure
  CV_size_group <- floor(n / CV_num_group)
  Omega_hat_CV_model_space <- array(0, dim=c(p, p, CV_num_group, num_model))
  Sigma_hat_CV_valid <- array(0, dim=c(p, p, CV_num_group))
  # The following loop takes a lot of time and memory.
  for (idx_CV_current_gourp in c(1:CV_num_group)) {
    # train for Omega, valid for Sigma
    X_train <- X[-c(((idx_CV_current_gourp-1)*CV_size_group+1):(idx_CV_current_gourp*CV_size_group)), ]
    X_valid <- X[c(((idx_CV_current_gourp-1)*CV_size_group+1):(idx_CV_current_gourp*CV_size_group)), ]
    # solution path is given
    Sigma_hat_CV_valid[ , , idx_CV_current_gourp] <- t(X_valid) %*% X_valid / dim(X_valid)[1]
    flag <- TRUE
    for (idx_curent_model in c(1:num_model)) {
      rs <- tryCatch(
        Omega_hat_CV_model_space[, , idx_CV_current_gourp, idx_curent_model] <- estimate_Omega(model_space[, idx_curent_model], t(X_train) %*% X_train / dim(X_train)[1]),
        error=function(e) {NULL}
      )
      if (is.null(rs)){
        tmp_num_model <- idx_curent_model - 1
        flag <- FALSE
      }
      if (!flag) break
    }
    if (!flag) num_model <- tmp_num_model # change loop condition variable
  }

  ## Estimate models using all samples
  Omega_hat_model_space <- array(0, dim=c(p, p, num_model))
  for (idx_curent_model in c(1:num_model)) {
    Omega_hat_model_space[, , idx_curent_model] <- estimate_Omega(model_space[, idx_curent_model], sample_covariance)
  }

  # # Choose Weight
  # Initial weight
  gamma_0 <- rep(1, num_model) / num_model

  CV_fun <- function(gamma_current) {
    # reture KL-loss under validation data
    res <- 0
    for (jj in c(1: CV_num_group)) {
      Omega_MA <- Reduce('+', lapply(c(1:num_model), function(ii) Omega_hat_CV_model_space[ , , jj, ii] * gamma_current[ii]))
      res <- res - log(det(Omega_MA)) + sum(diag(Omega_MA %*% Sigma_hat_CV_valid[ , , jj]))  # use covariance of validation samples or population_covariance
    }
    res
  }

  A <- rbind(diag(num_model), -diag(num_model))
  b <- matrix(c(rep(1, num_model), rep(0, num_model)), ncol=1)
  Aeq <- matrix(rep(1, num_model), nrow=1)
  beq <- 1
  gamma_hat <- pracma::fmincon(x0=gamma_0, fn=CV_fun, A=A, b=b, Aeq=Aeq, beq=beq)
  gamma_hat <- gamma_hat$par
  gamma_hat_plot <- gamma_hat

  Omega_hat_FMA <- Reduce('+', lapply(c(1:num_model), function(ii) Omega_hat_model_space[ , , ii] * gamma_hat[ii]))

  # output res: candidate models, weight, final estimator
  res <- list()
  res$model_set <- model_space
  res$weight <- gamma_hat
  res$Omega_hat <- Omega_hat_FMA

  return(res)
}