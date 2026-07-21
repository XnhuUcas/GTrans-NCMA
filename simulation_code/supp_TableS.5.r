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
p.list <- c(300) 
n.vec.list <- list(c(100, rep(225, K)))
Theta0.type='Random'


result_list <- list()
system.time({
  for (p in p.list) {
    for (n.vec in n.vec.list) {
      kl_loss_results <- list(
        kl_loss_GTrans_NCMA = numeric(n_repeats),
        kl_loss_trans_clime = numeric(n_repeats),
        kl_loss_clime = numeric(n_repeats),
        kl_loss_maegg = numeric(n_repeats),
        kl_loss_gnc_lasso = numeric(n_repeats),
        kl_loss_glasso = numeric(n_repeats),
        kl_loss_MT_glasso = numeric(n_repeats)
      )
      
      Frob_results <- list(
        Frob_GTrans_NCMA = numeric(n_repeats),
        Frob_trans_clime = numeric(n_repeats),
        Frob_clime = numeric(n_repeats),
        Frob_maegg = numeric(n_repeats),
        Frob_gnc_lasso = numeric(n_repeats),
        Frob_glasso = numeric(n_repeats),
        Frob_MT_glasso = numeric(n_repeats)
      )
      gamma_hat_results <- list(gamma_hat_GTrans_NCMA = list())
     
       time_results <- list(time_GTrans_NCMA = numeric(n_repeats),
                     time_trans_clime = numeric(n_repeats),
                     time_clime = numeric(n_repeats),
                     time_maegg = numeric(n_repeats),
                     time_gnc_lasso = numeric(n_repeats),
                     time_glasso = numeric(n_repeats),
                     time_MT_glasso = numeric(n_repeats)) 
      for (i in 1:n_repeats) {
        set.seed(100 + i)  # 设置不同的种子以产生不同的模拟数据
        
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
        
       ####### method 1: gnc_lasso_MA  (our method)######
        
        #generate Omega.list as candidate model for the final averaging process
        time_gncma <- system.time({
          Omega.list <- list()
          for (k in 1:length(new.X.list)) {
            
            X.std <- new.X.list[[k]]
            
            n <- nrow(X.std)
            
            Adj <- Data.list[[k]]$GridData$A
            d <- mean(colSums(Adj))
            
            alpha.cv<-sqrt(n)
            
            const<-0.5
            best.lambda <-const*2*sqrt(log(p)/nrow(X.std))
            
            t4.path <- GTrans_NCMA_single(
              X=X.std,
              A=Adj,
              lambda=best.lambda,
              alpha=alpha.cv
            )
            
            Omega_hat2 <- t4.path
            
            
            Omega.list[[k]] <- Omega_hat2
          }
          gncma_MA <- GTrans_NCMA(
            X=new.X.list[[1]],
            Adj=Data.list[[1]]$GridData$A,
            Omega.list=Omega.list,
            CV_num_group=5,
            V=5
          )
                })
        time_results$time_GTrans_NCMA[i] <- time_gncma[["elapsed"]]  
        Omega_hat_FMA<-gncma_MA$Omega_hat_FMA
        gamma_hat<-gncma_MA$gamma_hat

        res_kl2 <- -log(det(Omega_hat_FMA )) + sum(diag(Omega_hat_FMA  %*% Sigma_true)) - (-log(det(Omega_true)) + p)
        kl_loss_results$kl_loss_GTrans_NCMA[i] <- res_kl2
        Frob_results$Frob_GTrans_NCMA[i] <- sum((Omega_hat_FMA - Omega_true)^2) / p
        gamma_hat_results$gamma_hat_GTrans_NCMA[[i]] <- gamma_hat
        ####### method 2: trans_clime  ######
        time_trans_clime <- system.time({
        const=0.5
        n0=round(n.vec[1]*4/5) #split sample for aggregation
        Theta.re0<-Myfastclime.s(X=X.all[1:n0,], Bmat=diag(1,p), lambda=const*2*sqrt(log(p)/n0))
        Theta.init<-Theta.re0$Theta.hat
        Omega.tl1 <- Trans.CLIME(X=X.all[1:n0,], X.A=X.all[-(1:n.vec[1]),], const=const,
                                 X.til= X.all[(n0+1):n.vec[1],], Theta.cl=Theta.init)
        ind2<-(n.vec[1]-n0+1): n.vec[1]
        Theta.re0<-Myfastclime.s(X=X.all[ind2,], Bmat=diag(1,p), lambda=const*2*sqrt(log(p)/n0))
        Theta.init<-Theta.re0$Theta.hat
        Omega.tl2 <- Trans.CLIME(X=X.all[ind2,], X.A=X.all[-(1:n.vec[1]),], const=const,
                                 X.til= X.all[1:(n.vec[1]-n0),], Theta.cl=Theta.init)
        Omega.tl<-(Omega.tl1+Omega.tl2)/2
        })
        time_results$time_trans_clime[i] <- time_trans_clime[["elapsed"]]
        res_kl_trans_clime <- -log(det(Omega.tl)) + sum(diag(Omega.tl %*% Sigma_true)) - (-log(det(Omega_true)) + p)
        kl_loss_results$kl_loss_trans_clime[i] <- res_kl_trans_clime 
        Frob_results$Frob_trans_clime[i] <- sum((Omega.tl - Omega_true)^2) / p
        
        ####### method 3: MT_glasso  ######
        
        time_MT_glasso <- system.time({
        const.jgl<-jgl.fun(X.all,n.vec)$lam.const 
        jgl.re<-jgl.fun(X.all,n.vec, lam.const=const.jgl) 
        Omega.MT.glasso<-jgl.re$Theta.hat
        
        res_kl_MT.glasso <- -log(det(Omega.MT.glasso)) + sum(diag(Omega.MT.glasso %*% Sigma_true)) - (-log(det(Omega_true)) + p)
        })
        time_results$time_MT_glasso[i] <- time_MT_glasso[["elapsed"]]
        kl_loss_results$kl_loss_MT_glasso[i] <- res_kl_MT.glasso
        Frob_results$Frob_MT_glasso[i] <- sum((Omega.MT.glasso - Omega_true)^2) / p
        
        ####### method 4: maegg  ######
        
        # generate candidate models
        time_maegg <- system.time({
        glassopath_X <- glasso::glassopath(t(X.all[1:n.vec[1],])%*% (X.all[1:n.vec[1],])/(nrow(X.all[1:n.vec[1],])) , trace=FALSE)
        num_potential_model <- length(glassopath_X$rholist)
        # run maegg
        maegg.path <- maegg(X.all[1:n.vec[1],], CV_num_group=5, glassopath_X=glassopath_X, num_potential_model=num_potential_model, t(X.all[1:n.vec[1],])%*% (X.all[1:n.vec[1],])/(nrow(X.all[1:n.vec[1],])))
        Omega_hat_maegg <-  maegg.path$Omega_hat
        })
        time_results$time_maegg[i] <- time_maegg[["elapsed"]]
        res_kl_maegg <- -log(det(Omega_hat_maegg)) + sum(diag(Omega_hat_maegg %*% Sigma_true)) - (-log(det(Omega_true)) + p)
        kl_loss_results$kl_loss_maegg[i] <- res_kl_maegg
        Frob_results$Frob_maegg[i] <- sum((Omega_hat_maegg - Omega_true)^2) / p
        
        ####### method 5: gnc_lasso  ######
        time_gnc_lasso <- system.time({
        alpha.cv <-sqrt(n.vec[1])
        const<-0.5
        best.lambda <-const*2*sqrt(log(p)/n.vec[1])
        gnc_lasso.path <- gnclasso(X=X.all[1:n.vec[1],], A=Data.list[[1]]$GridData$A, lambda=best.lambda, alpha=alpha.cv)
        Omega_hat_gnc_lasso  <- solve(gnc_lasso.path$Sigma)
        })
        time_results$time_gnc_lasso[i] <- time_gnc_lasso[["elapsed"]]
        res_kl_gnc_lasso  <- -log(det(Omega_hat_gnc_lasso )) + sum(diag(Omega_hat_gnc_lasso  %*% Sigma_true)) - (-log(det(Omega_true)) + p)
        kl_loss_results$kl_loss_gnc_lasso[i] <- res_kl_gnc_lasso 
        Frob_results$Frob_gnc_lasso[i] <- sum((Omega_hat_gnc_lasso-Omega_true)^2) / p
        
        ####### method 6: clime ######
        time_clime <- system.time({
        Theta.re_clime<-Myfastclime.s(X=X.all[1:n.vec[1],], Bmat=diag(1,p), 
                                      lambda=const*2*sqrt(log(p)/n.vec[1]))
        Theta.hat_clime<-Theta.re_clime$Theta.hat
        })
        time_results$time_clime[i] <- time_clime[["elapsed"]]
        res_kl_clime<- -log(det(Theta.hat_clime)) + sum(diag(Theta.hat_clime %*% Sigma_true)) - (-log(det(Omega_true)) + p)
        kl_loss_results$kl_loss_clime[i] <- res_kl_clime
        Frob_results$Frob_clime[i] <- sum((Theta.hat_clime-Omega_true)^2) / p
        
        ####### method 7: glasso ######
        time_glasso <- system.time({
       X.sub <- X.all[1:n.vec[1], ]
X.sub <- scale(X.sub, center = TRUE, scale = FALSE)
S <- t(X.sub) %*% X.sub / nrow(X.sub)

lambda.g <- sqrt(log(p) / n.vec[1])
glasso.path <- glasso(s = S, rho = lambda.g, trace = FALSE)
          Omega_hat_glasso <- glasso.path$wi
        })
        time_results$time_glasso[i] <- time_glasso[["elapsed"]]  
        res_kl_glasso <- -log(det(Omega_hat_glasso)) + sum(diag(Omega_hat_glasso %*% Sigma_true)) - (-log(det(Omega_true)) + p)
        kl_loss_results$kl_loss_glasso[i] <- res_kl_glasso
        Frob_results$Frob_glasso[i] <- sum((Omega_hat_glasso - Omega_true)^2) / p
        
        
       
      }
      # Store the results for this specific p and n.vec combination
    result_list[[paste0("p_", p, "_nvec_", paste(n.vec, collapse = "_"))]] <- list(
  kl_loss = kl_loss_results,
  Frob = Frob_results,
  gamma_hat = gamma_hat_results,
  time = time_results
)
    }
  }
})

final_result<-result_list
