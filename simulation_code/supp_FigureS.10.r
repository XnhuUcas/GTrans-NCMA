library(glasso)
library(doParallel)
library(ggplot2)
library(reshape2)
library(RSpectra)
library(irlba)
library(lattice)
library(igraph)
library(clime)
library(CVglasso)
library(pracma)
library(NlcOptim)
library(lavaSearch2)
library(fastclime)
source("GNC-lasso.R")
source("maegg.R")
source("CLIME-opt.R")
source("method_functions.R")

n_repeats <-200
K=5
A.size=3
p.list <-  c(20,100) 
n.vec.list <- list(c(100, rep(225, K)),
                   c(225, rep(400, K)))
Theta0.type='Bdiag'


result_list <- list()
system.time({
  for (p in p.list) {
    for (n.vec in n.vec.list) {
      kl_loss_results <- list(
        kl_loss_maegg_pool     = numeric(n_repeats),
        kl_loss_gnc_lasso_pool = numeric(n_repeats),
        kl_loss_clime_pool     = numeric(n_repeats),
        kl_loss_glasso_pool    = numeric(n_repeats)
      )
      
      Frob_results <- list(
        Frob_maegg_pool     = numeric(n_repeats),
        Frob_gnc_lasso_pool = numeric(n_repeats),
        Frob_clime_pool     = numeric(n_repeats),
        Frob_glasso_pool    = numeric(n_repeats)
      )
      
      for (i in 1:n_repeats) {
        set.seed(100 + i) 
        
        # generate data
        sim.data<-sim.data <- DataGenerate(
          K=K,
          Info.size=A.size,
          h=0.4,
          n.vec=n.vec,
          s=10,
          p=p,
          type=Theta0.type)
        X <- sim.data$X
        Omega_true <- sim.data$Theta0
        Sigma_true <- solve(Omega_true)
        
        # split data into different categories
        X.split <- list()
        start_idx <- 1
        for (k in 1:(K + 1)) {
          end_idx <- start_idx + n.vec[k] - 1
          X.split[[k]] <- X[start_idx:end_idx, ]
          start_idx <- end_idx + 1
        }
        
        
        # generate mean matrices for every data set
        Data.list <- list()
        for (k in 1:length(X.split)) {
          
          X.class <- X.split[[k]]
          
          
          n <- nrow(X.class)
          
          m <- sqrt(n)
          GridData <- MeanmatrixGen(m)
          Adj <- GridData$A
          
          g <- graph.adjacency(Adj, mode = "undirected")
          D <- diag(colSums(Adj))
          L <- D - Adj
          
          M <- matrix(0, n, p)
          eig <- eigen(L, symmetric = TRUE)
          m1 <- eig$vectors[, n]
          m2 <- eig$vectors[, n - 1]
          mm <- sqrt(0.95) * m1 + sqrt(0.05) * m2
          
          for (j in 1:p) {
            M[, j] <- mm * sqrt(n)*1.5
          }
          Data.list[[k]] <- list(X.class = X.class, GridData = GridData, M = M)
        }
        
        #centralization for every data set
        new.X.list <- list()
        for (k in 1:length(X.split)) {
          X.class <- X.split[[k]]
          M <- Data.list[[k]]$M
          
          new.X <- X.class + M
          
          new.X.std <- scale(new.X, center = TRUE, scale = FALSE)
          
          new.X.list[[k]] <- new.X.std
        }
        
        X.all <- do.call(rbind,  new.X.list)
        ################################################################################################### 
        n_pool <- nrow(X.all)
        
        ############################
        ## Construct pooled A (block diagonal)
        ############################
        A_list <- lapply(Data.list, function(d) {
          as.matrix(d$GridData$A)
        })
        
        A_pool <- do.call(pracma::blkdiag, A_list)
        ################################################################################################### 
        
        ####################################################################
        ## Method 1: MAEGG-pool
        ####################################################################
        S_pool <- t(X.all) %*% X.all / n_pool
        
        glassopath_X <- glasso::glassopath(S_pool, trace = FALSE)
        num_potential_model <- length(glassopath_X$rholist)
        
        maegg.path <- maegg(
          X.all,
          CV_num_group = 5,
          glassopath_X = glassopath_X,
          num_potential_model = num_potential_model,
          S_pool
        )
        
        Omega_hat_maegg <- maegg.path$Omega_hat
        
        res_kl_maegg <- -log(det(Omega_hat_maegg)) +
          sum(diag(Omega_hat_maegg %*% Sigma_true)) -
          (-log(det(Omega_true)) + p)
        
        kl_loss_results$kl_loss_maegg_pool[i] <- res_kl_maegg
        Frob_results$Frob_maegg_pool[i] <- sum((Omega_hat_maegg - Omega_true)^2) / p
        
        ####################################################################
        ## Method 2: GNC-lasso-pool
        ####################################################################
        const <- 0.5
        alpha.cv <- sqrt(n_pool)
        best.lambda <- const * 2 * sqrt(log(p) / n_pool)
        
        gnc_lasso.path <- gnclasso(
          X = X.all,
          A = A_pool,
          lambda = best.lambda,
          alpha = alpha.cv
        )
        
        Omega_hat_gnc_lasso <- solve(gnc_lasso.path$Sigma)
        
        res_kl_gnc_lasso <- -log(det(Omega_hat_gnc_lasso)) +
          sum(diag(Omega_hat_gnc_lasso %*% Sigma_true)) -
          (-log(det(Omega_true)) + p)
        
        kl_loss_results$kl_loss_gnc_lasso_pool[i] <- res_kl_gnc_lasso
        Frob_results$Frob_gnc_lasso_pool[i] <- sum((Omega_hat_gnc_lasso - Omega_true)^2) / p
        
        ####################################################################
        ## Method 3: CLIME-pool
        ####################################################################
        Theta.re_clime <- Myfastclime.s(
          X = X.all,
          Bmat = diag(1, p),
          lambda = const * 2 * sqrt(log(p) / n_pool)
        )
        
        Theta.hat_clime <- Theta.re_clime$Theta.hat
        
        res_kl_clime <- -log(det(Theta.hat_clime)) +
          sum(diag(Theta.hat_clime %*% Sigma_true)) -
          (-log(det(Omega_true)) + p)
        
        kl_loss_results$kl_loss_clime_pool[i] <- res_kl_clime
        Frob_results$Frob_clime_pool[i] <- sum((Theta.hat_clime - Omega_true)^2) / p
        
        
        ####################################################################
        ## Method 4: glasso-pool
        ####################################################################
        glasso.path <- CVglasso(X.all)
        Omega_hat_glasso <- glasso.path$Omega
        
        res_kl_glasso <- -log(det(Omega_hat_glasso)) +
          sum(diag(Omega_hat_glasso %*% Sigma_true)) -
          (-log(det(Omega_true)) + p)
        
        kl_loss_results$kl_loss_glasso_pool[i] <- res_kl_glasso
        Frob_results$Frob_glasso_pool[i] <- sum((Omega_hat_glasso - Omega_true)^2) / p
        
        cat(
          " p =", p,
          "| n.vec =", paste(n.vec, collapse = ","),
          "| repeat =", i,
          "finished\n"
        )
        
      }
      # Store the results for this specific p and n.vec combination
      result_list[[paste0("p_", p, "_nvec_", paste(n.vec, collapse = "_"))]] <- list(
        kl_loss = kl_loss_results,
        Frob = Frob_results
      )
    }
  }
})

np_changing_result_list_type_Bdiag_pool <-result_list
