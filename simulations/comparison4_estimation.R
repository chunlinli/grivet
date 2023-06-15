# This file produces simulations corresponding to part (b) of Figure 1 in the supplementary materials.

## In the following `simu` functions, unless explicitly specified, 
## p: the number of primary variables.
## n: the number of observations.
## tau.list & gamma.list: lists of tunning parameters for TLP regression in parameter estimation.
## n.fold: number of folds for cross-validation in tuning parameters selection for causal effect matrix estimation.

library(doParallel)
n_cores <- detectCores() - 1
registerDoParallel(cores=n_cores)  
cl <- makeCluster(n_cores, type="FORK")

clusterEvalQ(cl,{library(grivet)
  library(mvtnorm)
  library(MASS)})

## This function apply our GrIVET to estimate the coefficients.

## Arguments:
## X: a n*q matrix for intervention variables.
## Y: a n*p matrix for primary variables.
## Pi: the matrix representing the support of ancestral relationships, 1 for existing ancestral relation, 0 otherwise.
## Phi: the matrix representing the support of interventional relationships, 1 for existing interventional relation, 0 otherwise.
## Piv: the matrix representing the candidate IV set, 1 for existing candidate IV relationship and 0 otherwise.
## tau.list & gamma.list: a list of tuning parameters to be selected for TLP regression.
## n.fold: the number of folds in selecting tuning parameters for TLP regression.

## Returns:
## U: the estimated causal effect matrix.
coef.proposed <- function(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold){
  out <- cv.intdag.coe(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  U <- out$U
  return(U)
}

## This function evaluate the performance of coefficient estimation.

## Arguments:
## U: the true causal effect matrix.
## U_es: the estimated causal effect matrix.
## Pi: the matrix representing true support of ancestral relationships, 1 for existing ancestral relation, 0 otherwise.

## Returns:
## max_abs_diff: maximum absolute difference between U and U_es.
## mean_abs_diff: mean absolute difference between U and U_es.
## mean_sq_diff: mean square difference between U and U_es.
metrics <- function(U,U_es,Pi){
  max_abs_diff <- max(abs(U-U_es))
  mean_abs_diff <- sum(abs(U-U_es))/max(1,sum(Pi!=0))
  mean_sq_diff <- sum((U-U_es)^2)/max(1,sum(Pi!=0))
  return(c(max_abs_diff, mean_abs_diff,mean_sq_diff))
}

## This function generate a hub graph with continuous intervention variables.

## Arguments: 
## p: the number of primary variables.
## q: the number of intervention variables.
## inv: the number of invalid instrumental variables for Y1.
## n: the number of observations.

## Returns:
## X: a n*q matrix for intervention variables.
## Y: a n*p matrix for primary variables.
## U: a p*P matrix representing the caual effect matrix.
## W: a q*p matrix representing the interventional effect matrix.
## Sigma: a p*p matrix representing the covariance matrix of residuals.

hub.generation <- function(p,q,inv,n){
  U <- matrix(0,p,p);U[1,2:p] <- 1;
  W1 <- diag(1,p,p); W2 <- matrix(1,inv,p); W3 <- matrix(0,q-p-inv,p); W3[,1] <- 1;W <- rbind(W1,W2,W3)
  vari <- diag(rep(0.5,p)); sigma <- vari%*%diag(p)%*%vari
  e <- mvrnorm(n,rep(0,p),sigma)
  r <- 1; V <- matrix(1,r,p); Z <- matrix(rnorm(r*n),ncol=r)
  X <- matrix(rnorm(q*n),ncol=q); Y <- (X%*%W+e+Z%*%V)%*%solve(diag(1,p)-U);
  Sigma <- sigma+t(V)%*%V
  return(list(X=X,Y=Y,U=U,W=W,Sigma=Sigma))
}

## This function generate a hub graph with discrete intervention variables.

## Arguments: 
## p: the number of primary variables.
## q: the number of intervention variables.
## inv: the number of invalid instrumental variables for Y1.
## n: the number of observations.

## Returns:
## X: a n*q matrix for intervention variables.
## Y: a n*p matrix for primary variables.
## U: a p*P matrix representing the caual effect matrix.
## W: a q*p matrix representing the interventional effect matrix.
## Sigma: a p*p matrix representing the covariance matrix of residuals.
hub.generation2 <- function(p,q,inv,n){
  U <- matrix(0,p,p);U[1,2:p] <- 1;
  W1 <- diag(1,p,p); W2 <- matrix(1,inv,p); W3 <- matrix(0,q-p-inv,p); W3[,1] <- 1;W <- rbind(W1,W2,W3)
  vari <- diag(rep(0.5,p)); sigma <- vari%*%diag(p)%*%vari
  e <- mvrnorm(n,rep(0,p),sigma)
  r <- 1; V <- matrix(1,r,p); Z <- matrix(rnorm(r*n),ncol=r)
  X <- matrix(sample(c(-1,1),q*n,replace = TRUE),ncol=q); Y <- (X%*%W+e+Z%*%V)%*%solve(diag(1,p)-U);
  Sigma <- sigma+t(V)%*%V
  return(list(X=X,Y=Y,U=U,W=W,Sigma=Sigma))
}

clusterExport(cl,c("coef.proposed","metrics","hub.generation","hub.generation2"))

## hub graph, n = 1000, continuous intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 10; inv <- 2; n <- 1000;
  tau.list <- seq(0.1,0.2,0.01)
  gamma.list <- seq(0.05,0.5,0.05)
  n.fold <- 5
  res <- numeric(18) # `res` is to store the metrics for the evaluation.
  
  ## Compute the evaluation metrics for GrIVET in parameter estimation with varying percentages of valid IVs. 
  ## u2: estimated coefficient matrix of DAG by GrIVET.
  for (q in 12:17){
    generation <- hub.generation(p,q,inv,n)
    X <- generation$X
    Y <- generation$Y
    U <- generation$U
    W <- generation$W
    V <- W%*%solve(diag(p)- U)
    out <- topological_order(V)
    Pi <- out$an_mat
    Phi <- out$in_mat
    Piv <- out$iv_mat
    u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
    res[((q-12)*3+1):((q-11)*3)] <- metrics(U,u2,Pi)
  }
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 18
clusterExport(cl,c("simu","len.stat"))
stats.hub <- parLapply(cl,1:1000,function(i){set.seed(i);stat <- simu(i);return(stat)})
stats_hub <- matrix(unlist(stats.hub),byrow = TRUE,ncol=len.stat)
write.csv(stats_hub,file.path("./primary_results/part4/estimation/","stats_hub1.csv"),row.names = FALSE)

## hub graph, n = 500, continuous intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 10; inv <- 2; n <- 500;
  tau.list <- seq(0.1,0.2,0.01)
  gamma.list <- seq(0.05,0.5,0.05)
  n.fold <- 5
  res <- numeric(18) # `res` is to store the metrics for the evaluation.
  
  ## Compute the evaluation metrics for GrIVET in parameter estimation with varying percentages of valid IVs. 
  ## u2: estimated coefficient matrix of DAG by GrIVET.
  for (q in 12:17){
    generation <- hub.generation(p,q,inv,n)
    X <- generation$X
    Y <- generation$Y
    U <- generation$U
    W <- generation$W
    V <- W%*%solve(diag(p)- U)
    out <- topological_order(V)
    Pi <- out$an_mat
    Phi <- out$in_mat
    Piv <- out$iv_mat
    u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
    res[((q-12)*3+1):((q-11)*3)] <- metrics(U,u2,Pi)
  }
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 18
clusterExport(cl,c("simu","len.stat"))
stats.hub <- parLapply(cl,1:1000,function(i){set.seed(i);stat <- simu(i);return(stat)})
stats_hub <- matrix(unlist(stats.hub),byrow = TRUE,ncol=len.stat)
write.csv(stats_hub,file.path("./primary_results/part4/estimation/","stats_hub2.csv"),row.names = FALSE)

## hub graph, n = 200, continuous intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 10; inv <- 2; n <- 200;
  tau.list <- seq(0.1,0.2,0.01)
  gamma.list <- seq(0.05,0.5,0.05)
  n.fold <- 5
  res <- numeric(18) # `res` is to store the metrics for the evaluation.
  
  ## Compute the evaluation metrics for GrIVET in parameter estimation with varying percentages of valid IVs. 
  ## u2: estimated coefficient matrix of DAG by GrIVET.
  for (q in 12:17){
    generation <- hub.generation(p,q,inv,n)
    X <- generation$X
    Y <- generation$Y
    U <- generation$U
    W <- generation$W
    V <- W%*%solve(diag(p)- U)
    out <- topological_order(V)
    Pi <- out$an_mat
    Phi <- out$in_mat
    Piv <- out$iv_mat
    u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
    res[((q-12)*3+1):((q-11)*3)] <- metrics(U,u2,Pi)
  }
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 18
clusterExport(cl,c("simu","len.stat"))
stats.hub <- parLapply(cl,1:1000,function(i){set.seed(i);stat <- simu(i);return(stat)})
stats_hub <- matrix(unlist(stats.hub),byrow = TRUE,ncol=len.stat)
write.csv(stats_hub,file.path("./primary_results/part4/estimation/","stats_hub3.csv"),row.names = FALSE)


## hub graph, n = 1000, discrete intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 10; inv <- 2; n <- 1000;
  tau.list <- seq(0.1,0.2,0.01)
  gamma.list <- seq(0.05,0.5,0.05)
  n.fold <- 5
  res <- numeric(18) # `res` is to store the metrics for the evaluation.
  
  ## Compute the evaluation metrics for GrIVET in parameter estimation with varying percentages of valid IVs. 
  ## u2: estimated coefficient matrix of DAG by GrIVET.
  for (q in 12:17){
    generation <- hub.generation2(p,q,inv,n)
    X <- generation$X
    Y <- generation$Y
    U <- generation$U
    W <- generation$W
    V <- W%*%solve(diag(p)- U)
    out <- topological_order(V)
    Pi <- out$an_mat
    Phi <- out$in_mat
    Piv <- out$iv_mat
    u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
    res[((q-12)*3+1):((q-11)*3)] <- metrics(U,u2,Pi)
  }
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 18
clusterExport(cl,c("simu","len.stat"))
stats.hub <- parLapply(cl,1:1000,function(i){set.seed(i);stat <- simu(i);return(stat)})
stats_hub <- matrix(unlist(stats.hub),byrow = TRUE,ncol=len.stat)
write.csv(stats_hub,file.path("./primary_results/part4/estimation/","stats_hub4.csv"),row.names = FALSE)

## hub graph, n = 500, discrete intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 10; inv <- 2; n <- 500;
  tau.list <- seq(0.1,0.2,0.01)
  gamma.list <- seq(0.05,0.5,0.05)
  n.fold <- 5
  res <- numeric(18) # `res` is to store the metrics for the evaluation.
  
  ## Compute the evaluation metrics for GrIVET in parameter estimation with varying percentages of valid IVs. 
  ## u2: estimated coefficient matrix of DAG by GrIVET.
  for (q in 12:17){
    generation <- hub.generation2(p,q,inv,n)
    X <- generation$X
    Y <- generation$Y
    U <- generation$U
    W <- generation$W
    V <- W%*%solve(diag(p)- U)
    out <- topological_order(V)
    Pi <- out$an_mat
    Phi <- out$in_mat
    Piv <- out$iv_mat
    u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
    res[((q-12)*3+1):((q-11)*3)] <- metrics(U,u2,Pi)
  }
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 18
clusterExport(cl,c("simu","len.stat"))
stats.hub <- parLapply(cl,1:1000,function(i){set.seed(i);stat <- simu(i);return(stat)})
stats_hub <- matrix(unlist(stats.hub),byrow = TRUE,ncol=len.stat)
write.csv(stats_hub,file.path("./primary_results/part4/estimation/","stats_hub5.csv"),row.names = FALSE)

## hub graph, n = 200, discrete intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 10; inv <- 2; n <- 200;
  tau.list <- seq(0.1,0.2,0.01)
  gamma.list <- seq(0.05,0.5,0.05)
  n.fold <- 5
  res <- numeric(18) # `res` is to store the metrics for the evaluation.
  
  ## Compute the evaluation metrics for GrIVET in parameter estimation with varying percentages of valid IVs. 
  ## u2: estimated coefficient matrix of DAG by GrIVET.
  for (q in 12:17){
    generation <- hub.generation2(p,q,inv,n)
    X <- generation$X
    Y <- generation$Y
    U <- generation$U
    W <- generation$W
    V <- W%*%solve(diag(p)- U)
    out <- topological_order(V)
    Pi <- out$an_mat
    Phi <- out$in_mat
    Piv <- out$iv_mat
    u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
    res[((q-12)*3+1):((q-11)*3)] <- metrics(U,u2,Pi)
  }
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 18
clusterExport(cl,c("simu","len.stat"))
stats.hub <- parLapply(cl,1:1000,function(i){set.seed(i);stat <- simu(i);return(stat)})
stats_hub <- matrix(unlist(stats.hub),byrow = TRUE,ncol=len.stat)
write.csv(stats_hub,file.path("./primary_results/part4/estimation/","stats_hub6.csv"),row.names = FALSE)


stopCluster(cl)
