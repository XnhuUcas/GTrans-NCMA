library(glasso)
library(mvtnorm)


MeanmatrixGen <- function(m){
  n <- m^2
  Adj <- matrix(0, n, n)
  value <- rep(0, n)
  x <- 1:m
  value.mat <- outer(x, x, function(a, b) (a / m) * (b / m))

  for(i in 1:m){
    for(j in 1:m){
      v <- (j - 1) * m + i  
      value[v] <- value.mat[i, j]
      if(i > 1){
        v.upper <- (j - 1) * m + (i - 1)
        Adj[v, v.upper] <- 1
      }
      if(i < m){
        v.lower <- (j - 1) * m + (i + 1)
        Adj[v, v.lower] <- 1
      }
      if(j > 1){
        v.left <- (j - 2) * m + i
        Adj[v, v.left] <- 1
      }
      if(j < m){
        v.right <- j * m + i
        Adj[v, v.right] <- 1
      }
    }
  }

  return(list(A = Adj, value = value, value.mat = value.mat))
}
                     

GTrans_NCMA_B_semi <- function(X,A,alpha,L.U=NULL,L.tau=NULL,low.rank=NULL,train.index){
  A <- as.matrix(A)
  D <- diag(colSums(A))
  L <- D-A
  n <- nrow(A)
  p <- ncol(X)
  if(ncol(A)!=n) stop("Invalid adjacency matrix!")
  
  if(is.null(L.U)){
    if(is.null(low.rank)){
      L.eigen <- eigen(L)
      L.U <- L.eigen$vectors
      L.tau <- L.eigen$values + 0.01*min(L.eigen$values)
    }else{
      L.eigen <-	partial_eigen(L, n = low.rank, symmetric = TRUE)
      L.tau <- L.eigen$values + 0.01*min(L.eigen$values)
      L.U <- L.eigen$vectors
    }
  }
  X.train <- X[train.index,]
  X.test <- X[-train.index,]
  test.index <- setdiff(1:nrow(X),train.index)
  L.U.train <- L.U[train.index,]
  M <- solve(t(L.U.train)%*%L.U.train+alpha*diag(L.tau),t(L.U.train)%*%X.train)
  B <- L.U.train%*%M
  B.test <- L.U[test.index,]%*%M
  B.all <- matrix(0,nrow(X),p)
  B.all[train.index,] <- B
  B.all[test.index,] <- B.test
  
  err <- norm(B.test-X.test,"F")^2/(length(test.index)*p)
  return(list(B.train=B,B.test=B.test,err=err,B.all=B.all))
}


DataGenerate <- function(K=10, Info.size=4, h, n.vec, s=10, p=100, type='Toep') {
  if(type=='Toep'){
    Theta <- toeplitz(0.6^(1:p)*2)
    Theta[which(abs(Theta)<=0.05, arr.ind=T)] <- 0
  } else if(type=='Bdiag'){
    Theta <- kronecker(diag(p/4), toeplitz(c(1.2,0.9,0.6,0.3)))
  } else if(type == 'Random'){
    Theta <- diag(1,p) + matrix(runif(p^2,0, 0.8),ncol=p,nrow=p)
    Theta <- (Theta + t(Theta))/2
    for(j in 1:p){
      for(l in 1:p){
        Theta[j,l] = Theta[j,l]/sqrt(abs(j-l)+1)
      }
    }
    for(j in 1:p){
      Theta[,j] <- Theta[,j]*(abs(Theta[,j])>=quantile(abs(Theta[,j]),1-s/p))
      Theta[j,] <- Theta[j,]*(abs(Theta[j,])>=quantile(abs(Theta[j,]),1-s/p))
    }
    Theta <- diag(max(0.1-min(eigen(Theta)$values),0),p) + Theta
  }
  Sig_correct <- solve(Theta)
  if(min(eigen(Sig_correct)$values)<0.05){
    Sig_correct <- Sig_correct + diag(0.1-min(eigen(Sig_correct)$values),p)
  }
  Sig <- solve(Theta)
  X <- rmvnorm(n.vec[1], rep(0,p), sigma=Sig_correct)
  Theta.out <- diag(1,p)
  Omega.vec <- 0
  Sig_inv_list_source <- list()
  for(k in 1 : K){
    if(k < Info.size+1){
      Delta.k <- matrix(rbinom(p^2,size=1,prob=0.1)*runif(p^2,-h/p,h/p),ncol=p)
      Sig.k <- (diag(1,p)+Delta.k)%*%Sig
      Sig.k <- (Sig.k+t(Sig.k))/2
      if(min(eigen(Sig.k)$values)<0.05){
        Sig.k <- Sig.k+diag(0.1-min(eigen(Sig.k)$values),p)
      }
      X <- rbind(X, rmvnorm(n.vec[k+1],rep(0,p), sigma=Sig.k))
      Sig_inv_list_source[[k]] <- solve(Sig.k)
    }
    if(Info.size < k && k < Info.size+3){
      Delta.out <- diag(1,p)+matrix(rbinom(p^2,size=1,prob=0.3)*0.7,ncol=p)
      Sig.k <- (diag(1,p)+Delta.out)%*%Sig
      Sig.k <- (Sig.k+t(Sig.k))/2
      if(min(eigen(Sig.k)$values)<0.05){
        Sig.k <- Sig.k+diag(0.1-min(eigen(Sig.k)$values),p)
      }
      X <- rbind(X, rmvnorm(n.vec[k+1],rep(0,p), sigma=Sig.k))
      Sig_inv_list_source[[k]] <- solve(Sig.k)
    }
    Omega.vec <- c(Omega.vec, max(colSums(abs(diag(1,p)-Sig.k%*%Theta))))
  }
  cat(Omega.vec,'\n')
  list(X=X, Theta0=Theta, Omega.l1=max(Omega.vec))
}



estimate_Omega_constrained <- function(graphical_model, covariance_train){
  p <- dim(covariance_train)[1]
  Omega_hat <- matrix(rep(0, p*p), p)
  beta_hat <- matrix(rep(0, p*p),p)
  adjacency_matrix <- matrix(rep(0, p*p), p)
  adjacency_matrix[upper.tri(adjacency_matrix)] <- graphical_model
  adjacency_matrix <- adjacency_matrix + t(adjacency_matrix) + diag(p)
  convergence_status <- rep(FALSE, p)
  current_idx <- 0
  W_current <- W_next <- covariance_train
  iter_num <- 0
  while (!identical(convergence_status, rep(TRUE, p))) {
    current_idx <- current_idx + 1 - (current_idx %/% p)*p
    adjacency_nodes_with_current <- setdiff(which(adjacency_matrix[current_idx,] != 0), current_idx)
    W11 <- W_next[-current_idx, -current_idx]
    W11_star <- W_next[adjacency_nodes_with_current, adjacency_nodes_with_current]
    s12_star <- covariance_train[current_idx, adjacency_nodes_with_current]
    if (length(adjacency_nodes_with_current) > 0){
      beta_star <- solve(as.matrix(W11_star)) %*% s12_star
      beta_hat[adjacency_nodes_with_current, current_idx] <- beta_star
    }
    W12 <- W11 %*% beta_hat[-current_idx, current_idx]
    W_next[-current_idx, current_idx] <- W12
    iter_num <- iter_num + 1
    if (norm(W_current-W_next, "F") < 1e-5) {
      convergence_status[current_idx] = TRUE
      W_hat <- W_next
    }
    W_current <- W_next
  }
  for (current_idx in c(1:p)) {
    s22 <- covariance_train[current_idx, current_idx]
    W12 <- W_hat[-current_idx, current_idx]
    Omega22 <- 1 / (s22 - W12 %*% beta_hat[-current_idx, current_idx])
    Omega12 <- - beta_hat[-current_idx,current_idx] * as.numeric(Omega22)
    Omega_hat[-current_idx, current_idx] <- Omega12
    Omega_hat[current_idx, -current_idx] <- Omega12
    Omega_hat[current_idx, current_idx] <- Omega22
  }
  gc()
  return(Omega_hat)
}

GTrans_NCMA_B_cv <- function(X,A,alpha.seq,V,proportion=0.1,L.U=NULL,L.tau=NULL,low.rank=NULL){
  n <- nrow(A)
  p <- ncol(X)
  K <- V
  m <- floor(nrow(X)*proportion)
  A <- as.matrix(A)
  D <- diag(colSums(A))
  L <- D-A
  if(ncol(A)!=n) stop("Invalid adjacency matrix!")
  if(is.null(L.U)){
    if(is.null(low.rank)){
      L.eigen <- eigen(L)
      L.U <- L.eigen$vectors
      L.tau <- L.eigen$values + 0.00001*min(L.eigen$values)
    }else{
      L.eigen <- partial_eigen(L, n = low.rank, symmetric = TRUE)
      L.tau <- L.eigen$values + 0.00001*min(L.eigen$values)
      L.U <- L.eigen$vectors
    }
  }
  cv.err.mat <- matrix(0,K,length(alpha.seq))
  for(k in 1:K){
    print(paste("CV iteration: ",k))
    train.index <- sample(nrow(X),m)
    for(i in 1:length(alpha.seq)){
      cv.semi <- GTrans_NCMA_B_semi(X,A,alpha.seq[i],L.U=L.U,L.tau=L.tau,low.rank=low.rank,train.index=train.index)
      cv.err.mat[k,i] <- cv.semi$err
    }
  }
  return(list(cv.err.mat=cv.err.mat,opt.index=which.min(colMeans(cv.err.mat))))
}

GTrans_NCMA_B <- function(X,A,alpha,X_index,L.U=NULL,L.tau=NULL,low.rank=NULL){
  A <- as.matrix(A)
  D <- diag(colSums(A))
  L <- D-A
  n <- nrow(A)
  p <- ncol(X)
  if(ncol(A)!=n) stop("Invalid adjacency matrix!")
  if(is.null(L.U)){
    if(is.null(low.rank)){
      L.eigen <- eigen(L)
      L.U <- L.eigen$vectors
      L.tau <- L.eigen$values + 0.00001*min(L.eigen$values)
    }else{
      L.eigen <- partial_eigen(L, n = low.rank, symmetric = TRUE)
      L.tau <- L.eigen$values + 0.00001*min(L.eigen$values)
      L.U <- L.eigen$vectors
    }
  }
  B.init <- matrix(0,nrow=nrow(X),ncol=p)
  L.U.train <- L.U[X_index,]
  M <- solve(t(L.U.train)%*%L.U.train+alpha*diag(L.tau),t(L.U.train)%*%X)
  B.init <- L.U.train%*%M
  B <- B.init
  return(B)
}

GTrans_NCMA_single <- function(X, A, lambda= NULL, alpha = 1,tol=1e-4,verbose=FALSE,W.init=NULL,Theta.init=NULL,M.init=NULL,L.U=NULL,L.tau=NULL,low.rank=NULL) {
  p <- ncol(X)
  n <- nrow(X)
  B_hat_all <- GTrans_NCMA_B(X,A, alpha ,L.U=NULL,L.tau=NULL,low.rank=NULL)
  S <- t(X-B_hat_all)%*%(X-B_hat_all)/n
  glasso_results <- glasso(S, rho = lambda, trace = FALSE)
  Theta_glasso <- glasso_results$wi
  support <- (abs(Theta_glasso) > 1e-8) * 1
  diag(support) <- 1
  support <- support[upper.tri(support)]
  Omega_hat <- estimate_Omega_constrained(support, S)
  return(
    Omega = Omega_hat
  )
}

GTrans_NCMA_single_cv <- function(X, X_index, A, lambda, alpha,
                                  tol = 1e-4, verbose = FALSE,
                                  W.init = NULL, Theta.init = NULL,
                                  M.init = NULL, L.U = NULL,
                                  L.tau = NULL, low.rank = NULL) {
  A <- as.matrix(A)
  D <- diag(colSums(A))
  L <- D-A
  n <- nrow(A)
  p <- ncol(X)
  if(ncol(A)!=n) stop("Invalid adjacency matrix!")
  if(is.null(W.init)){
    S <- t(X)%*%X/n
    diag(S) <- diag(S)+lambda
    W.init <- S
  }
  if(is.null(Theta.init)){
    Theta.init <- solve(W.init)
  }
  if(is.null(L.U)){
    if(is.null(low.rank)){
      L.eigen <- eigen(L)
      L.U <- L.eigen$vectors
      L.tau <- L.eigen$values + 0.01*min(L.eigen$values)
    }else{
      L.eigen <- partial_eigen(L, n = low.rank, symmetric = TRUE)
      L.tau <- L.eigen$values + 0.01*min(L.eigen$values)
      L.U <- L.eigen$vectors
    }
  }
  if(is.null(M.init)){
    M.init <- matrix(0,nrow=n,ncol=p)
    L.U.train <- L.U[X_index,]
    M <- solve(t(L.U.train)%*%L.U.train+alpha*diag(L.tau),t(L.U.train)%*%X)
    B.init <- L.U.train%*%M
  }
  B <- B.init
  B.centered <- X-B
  S <- t(B.centered)%*%B.centered/n
  glasso_results <- glasso(S, rho = lambda, trace = FALSE)
  Theta_glasso <- glasso_results$wi
  support <- (abs(Theta_glasso) > 1e-8) * 1
  diag(support) <- 1
  support <- support[upper.tri(support)]
  result.g <- estimate_Omega_constrained(support, S)
  result <- list()
  result$B <- B
  result$Sigma <- result.g
  return(result)
}

GTrans_NCMA <- function(X, Adj, Omega.list, CV_num_group, V, alpha.seq = NULL) {
  p <- ncol(X)
  n <- nrow(X)
  K <- length(Omega.list) - 1
  Omega_hat_CV_model_space <- array(0, dim = c(p, p, CV_num_group, K + 1))
  Sigma_hat_CV_valid <- array(0, dim = c(p, p, CV_num_group))
  CV_size_group <- floor(n / CV_num_group)
  for (idx_CV_current_group in 1:CV_num_group) {
    X_train_first <- X[-c(((idx_CV_current_group - 1) * CV_size_group + 1):(idx_CV_current_group * CV_size_group)), ]
    X_valid_first <- X[c(((idx_CV_current_group - 1) * CV_size_group + 1):(idx_CV_current_group * CV_size_group)), ]
    Omega_train_rest <- Omega.list[2:(K + 1)]
    X_train_first_std <- scale(X_train_first, center=TRUE, scale=FALSE)
    X_valid_first_std <- scale(X_valid_first, center=TRUE, scale=FALSE)
    rho.seq <- exp(seq(log(0.002), log(1), length.out = 20))
    alpha.seq <- rho.seq * n / d
    test.B.cv <- GTrans_NCMA_B_cv(X_train_first_std, Adj, alpha.seq, V = V)
    alpha.cv <- alpha.seq[test.B.cv$opt.index]
    const <- 0.5
    best.lambda <- const*2*sqrt(log(p)/nrow(X_train_first_std))
    t90.path <- GTrans_NCMA_single_cv(X_train_first_std,-c(((idx_CV_current_group-1)*CV_size_group+1):(idx_CV_current_group*CV_size_group)), A= Adj, lambda= best.lambda, alpha= alpha.cv)
    Omega_train_first <- t90.path$Sigma
    Omega_hat_CV_model_space[ , , idx_CV_current_group, 1] <- Omega_train_first
    for (k in 2:(K+1)) {
      Omega_hat_CV_model_space[ , , idx_CV_current_group, k] <- Omega_train_rest[[k-1]]
    }
    cv_valid <- GTrans_NCMA_B_cv(X_valid_first, Adj, alpha.seq, V = V, proportion = 0.2)
    cv_valid_opt.index <- alpha.seq[cv_valid$opt.index]
    B_hat_valid <- GTrans_NCMA_B(X_valid_first, Adj, cv_valid_opt.index, c(((idx_CV_current_group-1)*CV_size_group+1):(idx_CV_current_group*CV_size_group)))
    Sigma_hat_CV_valid[ , , idx_CV_current_group] <- t(X_valid_first - B_hat_valid) %*% (X_valid_first - B_hat_valid) / nrow(X_valid_first)
  }
  gamma_0 <- rep(1, K + 1) / (K + 1)
  CV_fun <- function(gamma) {
    sum(sapply(1:CV_num_group, function(jj) {
      Omega_MA <- Reduce("+", lapply(1:(K+1), function(ii) {
        Omega_hat_CV_model_space[, , jj, ii] * gamma[ii]
      }))
      -log(det(Omega_MA)) + sum(diag(Omega_MA %*% Sigma_hat_CV_valid[, , jj]))
    }))
  }
  Q <- rbind(diag(K+1), -diag(K+1))
  b <- matrix(c(rep(1, K+1), rep(0, K+1)), ncol = 1)
  Aeq <- matrix(rep(1, K+1), nrow = 1)
  beq <- 1
  gamma_hat <- pracma::fmincon(gamma_0, CV_fun, A = Q, b = b, Aeq = Aeq, beq = beq)$par
  Omega_hat_FMA <- Reduce("+", lapply(1:(K+1), function(ii) {
    Omega.list[[ii]] * gamma_hat[ii]
  }))
  return(list(
    Omega_hat = Omega_hat_FMA,
    gamma_hat = gamma_hat,
    tau = sum(gamma_hat[c(1, 3:7, 9:10)])
  ))
}