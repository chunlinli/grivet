# This file produces simulations corresponding to Table 2 in the main text.

## In the following `simu` functions, unless explicitly specified, 
## p: the number of primary variables.
## n: the number of observations.
## tau.list & gamma.list: lists of tunning parameters for TLP regression in parameter estimation.
## n.fold: number of folds for cross-validation in tuning parameters selection for causal effect matrix estimation.
## type: type = 0 controls the simulated graph so that H_0 holds, otherwise H_1 holds.
library(doParallel)
n_cores <- detectCores() - 1
registerDoParallel(cores=n_cores)  
cl <- makeCluster(n_cores, type="FORK")

clusterEvalQ(cl,{source("comparison2_setting1.R",chdir = TRUE)
  library(mvtnorm)
  library(MASS)})

## random graph, n = 500, continuous intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 100;n <- 500;
  tau.list <- seq(0.2,0.3,0.01)
  gamma.list <- seq(0.1,1,0.1)
  n.fold <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(6) # `res` is to store the metrics for the evaluation.
  
  ## Data generation, Y stands for the primary variables and X stands for the intervention variables.
  ## Refer to descriptions of simulations related to Table 2 in the main text for details.
  generation <- random.generation(p,n,U.test,type)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  W <- generation$W
  V <- W%*%solve(diag(p)- U)
  
  
  out <- topological_order(V)
  Pi <- out$an_mat ## The matrix representing ancestral relationship, 1 if existing, 0 otherwise.
  Phi <- out$in_mat ## The matrix representing invertion relationship, 1 if existing, 0 otherwise.
  Piv <- out$iv_mat ## The matrix representing candidate IV set, 1 if existing, 0 otherwise.
  
  
  ## Compute the evaluation metrics for method without taking consideration of confounders in parameter estimation. 
  ## u1: estimated coefficient matrix of DAG by directly sparse regression.
  u1 <- coef.direct(X,Y,Pi,Phi,tau.list,gamma.list,n.fold)
  res[1:3] <- metrics(U,u1,Pi)
  
  ## Compute the evaluation metrics for GrIVET in parameter estimation. 
  ## u2: estimated coefficient matrix of DAG by GrIVET.
  u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  res[4:6] <- metrics(U,u2,Pi)
  
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 6
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./primary_results/part2/setting1/","stats_random1.csv"),row.names = FALSE)

## random graph, n = 400, continuous intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 100;n <- 400;
  tau.list <- seq(0.2,0.3,0.01)
  gamma.list <- seq(0.1,1,0.1)
  n.fold <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(6) # `res` is to store the metrics for the evaluation.
  
  ## Data generation, Y stands for the primary variables and X stands for the intervention variables.
  ## Refer to descriptions of simulations related to Table 2 in the main text for details.
  generation <- random.generation(p,n,U.test,type)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  W <- generation$W
  V <- W%*%solve(diag(p)- U)
  
  out <- topological_order(V)
  Pi <- out$an_mat ## The matrix representing ancestral relationship, 1 if existing, 0 otherwise.
  Phi <- out$in_mat ## The matrix representing invertion relationship, 1 if existing, 0 otherwise.
  Piv <- out$iv_mat ## The matrix representing candidate IV set, 1 if existing, 0 otherwise.
  
  
  ## Compute the evaluation metrics for method without taking consideration of confounders in parameter estimation. 
  ## u1: estimated coefficient matrix of DAG by directly sparse regression.
  u1 <- coef.direct(X,Y,Pi,Phi,tau.list,gamma.list,n.fold)
  res[1:3] <- metrics(U,u1,Pi)
  
  ## Compute the evaluation metrics for GrIVET in parameter estimation. 
  ## u2: estimated coefficient matrix of DAG by GrIVET.
  u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  res[4:6] <- metrics(U,u2,Pi)
  
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 6
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./primary_results/part2/setting1/","stats_random2.csv"),row.names = FALSE)

## random graph, n = 300, continuous intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 100;n <- 300;
  tau.list <- seq(0.2,0.3,0.01)
  gamma.list <- seq(0.1,1,0.1)
  n.fold <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(6) # `res` is to store the metrics for the evaluation.
  
  ## Data generation, Y stands for the primary variables and X stands for the intervention variables.
  ## Refer to descriptions of simulations related to Table 2 in the main text for details.
  generation <- random.generation(p,n,U.test,type)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  W <- generation$W
  V <- W%*%solve(diag(p)- U)
  
  out <- topological_order(V)
  Pi <- out$an_mat ## The matrix representing ancestral relationship, 1 if existing, 0 otherwise.
  Phi <- out$in_mat ## The matrix representing invertion relationship, 1 if existing, 0 otherwise.
  Piv <- out$iv_mat ## The matrix representing candidate IV set, 1 if existing, 0 otherwise.
  
  
  ## Compute the evaluation metrics for method without taking consideration of confounders in parameter estimation. 
  ## u1: estimated coefficient matrix of DAG by directly sparse regression.
  u1 <- coef.direct(X,Y,Pi,Phi,tau.list,gamma.list,n.fold)
  res[1:3] <- metrics(U,u1,Pi)
  
  ## Compute the evaluation metrics for GrIVET in parameter estimation. 
  ## u2: estimated coefficient matrix of DAG by GrIVET.
  u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  res[4:6] <- metrics(U,u2,Pi)
  
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 6
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./primary_results/part2/setting1/","stats_random3.csv"),row.names = FALSE)

## random graph, n = 500, discrete intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 100;n <- 500;
  tau.list <- seq(0.2,0.3,0.01)
  gamma.list <- seq(0.1,1,0.1)
  n.fold <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(6) # `res` is to store the metrics for the evaluation.
  
  ## Data generation, Y stands for the primary variables and X stands for the intervention variables.
  ## Refer to descriptions of simulations related to Table 2 in the main text for details.
  generation <- random.generation2(p,n,U.test,type)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  W <- generation$W
  V <- W%*%solve(diag(p)- U)
  
  out <- topological_order(V)
  Pi <- out$an_mat ## The matrix representing ancestral relationship, 1 if existing, 0 otherwise.
  Phi <- out$in_mat ## The matrix representing invertion relationship, 1 if existing, 0 otherwise.
  Piv <- out$iv_mat ## The matrix representing candidate IV set, 1 if existing, 0 otherwise.
  
  
  ## Compute the evaluation metrics for method without taking consideration of confounders in parameter estimation. 
  ## u1: estimated coefficient matrix of DAG by directly sparse regression.
  u1 <- coef.direct(X,Y,Pi,Phi,tau.list,gamma.list,n.fold)
  res[1:3] <- metrics(U,u1,Pi)
  
  ## Compute the evaluation metrics for GrIVET in parameter estimation. 
  ## u2: estimated coefficient matrix of DAG by GrIVET.
  u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  res[4:6] <- metrics(U,u2,Pi)
  
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 6
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./primary_results/part2/setting1/","stats_random4.csv"),row.names = FALSE)

## random graph, n = 400, discrete intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 100;n <- 400;
  tau.list <- seq(0.2,0.3,0.01)
  gamma.list <- seq(0.1,1,0.1)
  n.fold <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(6) # `res` is to store the metrics for the evaluation.
  
  ## Data generation, Y stands for the primary variables and X stands for the intervention variables.
  ## Refer to descriptions of simulations related to Table 2 in the main text for details.
  generation <- random.generation2(p,n,U.test,type)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  W <- generation$W
  V <- W%*%solve(diag(p)- U)
  
  out <- topological_order(V)
  Pi <- out$an_mat ## The matrix representing ancestral relationship, 1 if existing, 0 otherwise.
  Phi <- out$in_mat ## The matrix representing invertion relationship, 1 if existing, 0 otherwise.
  Piv <- out$iv_mat ## The matrix representing candidate IV set, 1 if existing, 0 otherwise.
  
  
  ## Compute the evaluation metrics for method without taking consideration of confounders in parameter estimation. 
  ## u1: estimated coefficient matrix of DAG by directly sparse regression.
  u1 <- coef.direct(X,Y,Pi,Phi,tau.list,gamma.list,n.fold)
  res[1:3] <- metrics(U,u1,Pi)
  
  ## Compute the evaluation metrics for GrIVET in parameter estimation. 
  ## u2: estimated coefficient matrix of DAG by GrIVET.
  u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  res[4:6] <- metrics(U,u2,Pi)
  
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 6
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./primary_results/part2/setting1/","stats_random5.csv"),row.names = FALSE)

## random graph, n = 300, discrete intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 100;n <- 300;
  tau.list <- seq(0.2,0.3,0.01)
  gamma.list <- seq(0.1,1,0.1)
  n.fold <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(6) # `res` is to store the metrics for the evaluation.
  
  ## Data generation, Y stands for the primary variables and X stands for the intervention variables.
  ## Refer to descriptions of simulations related to Table 2 in the main text for details.
  generation <- random.generation2(p,n,U.test,type)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  W <- generation$W
  V <- W%*%solve(diag(p)- U)
  
  out <- topological_order(V)
  Pi <- out$an_mat ## The matrix representing ancestral relationship, 1 if existing, 0 otherwise.
  Phi <- out$in_mat ## The matrix representing invertion relationship, 1 if existing, 0 otherwise.
  Piv <- out$iv_mat ## The matrix representing candidate IV set, 1 if existing, 0 otherwise.
  
  
  ## Compute the evaluation metrics for method without taking consideration of confounders in parameter estimation. 
  ## u1: estimated coefficient matrix of DAG by directly sparse regression.
  u1 <- coef.direct(X,Y,Pi,Phi,tau.list,gamma.list,n.fold)
  res[1:3] <- metrics(U,u1,Pi)
  
  ## Compute the evaluation metrics for GrIVET in parameter estimation. 
  ## u2: estimated coefficient matrix of DAG by GrIVET.
  u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  res[4:6] <- metrics(U,u2,Pi)
  
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 6
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./primary_results/part2/setting1/","stats_random6.csv"),row.names = FALSE)


## hub graph, n = 500, continuous intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 101;n <- 500;
  tau.list <- seq(0.2,0.3,0.01)
  gamma.list <- seq(0.1,1,0.1)
  n.fold <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(6) # `res` is to store the metrics for the evaluation.
  
  ## Data generation, Y stands for the primary variables and X stands for the intervention variables.
  ## Refer to descriptions of simulations related to Table 2 in the main text for details.
  generation <- hub.generation(p,n)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  W <- generation$W
  V <- W%*%solve(diag(p)- U)
  
  out <- topological_order(V)
  Pi <- out$an_mat ## The matrix representing ancestral relationship, 1 if existing, 0 otherwise.
  Phi <- out$in_mat ## The matrix representing invertion relationship, 1 if existing, 0 otherwise.
  Piv <- out$iv_mat ## The matrix representing candidate IV set, 1 if existing, 0 otherwise.
  
  
  ## Compute the evaluation metrics for method without taking consideration of confounders in parameter estimation. 
  ## u1: estimated coefficient matrix of DAG by directly sparse regression.
  u1 <- coef.direct(X,Y,Pi,Phi,tau.list,gamma.list,n.fold)
  res[1:3] <- metrics(U,u1,Pi)
  
  ## Compute the evaluation metrics for GrIVET in parameter estimation. 
  ## u2: estimated coefficient matrix of DAG by GrIVET.
  u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  res[4:6] <- metrics(U,u2,Pi)
  
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 6
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);stat <- simu(i);return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./primary_results/part2/setting1/","stats_hub1.csv"),row.names = FALSE)

## hub graph, n = 400, continuous intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 101;n <- 400;
  tau.list <- seq(0.2,0.3,0.01)
  gamma.list <- seq(0.1,1,0.1)
  n.fold <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(6) # `res` is to store the metrics for the evaluation.
  
  ## Data generation, Y stands for the primary variables and X stands for the intervention variables.
  ## Refer to descriptions of simulations related to Table 2 in the main text for details.
  generation <- hub.generation(p,n)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  W <- generation$W
  V <- W%*%solve(diag(p)- U)
  
  out <- topological_order(V)
  Pi <- out$an_mat ## The matrix representing ancestral relationship, 1 if existing, 0 otherwise.
  Phi <- out$in_mat ## The matrix representing invertion relationship, 1 if existing, 0 otherwise.
  Piv <- out$iv_mat ## The matrix representing candidate IV set, 1 if existing, 0 otherwise.
  
  
  ## Compute the evaluation metrics for method without taking consideration of confounders in parameter estimation. 
  ## u1: estimated coefficient matrix of DAG by directly sparse regression.
  u1 <- coef.direct(X,Y,Pi,Phi,tau.list,gamma.list,n.fold)
  res[1:3] <- metrics(U,u1,Pi)
  
  ## Compute the evaluation metrics for GrIVET in parameter estimation. 
  ## u2: estimated coefficient matrix of DAG by GrIVET.
  u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  res[4:6] <- metrics(U,u2,Pi)
  
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 6
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);stat <- simu(i);return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./primary_results/part2/setting1/","stats_hub2.csv"),row.names = FALSE)


## hub graph, n = 300, continuous intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 101;n <- 300;
  tau.list <- seq(0.2,0.3,0.01)
  gamma.list <- seq(0.1,1,0.1)
  n.fold <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(6) # `res` is to store the metrics for the evaluation.
  
  ## Data generation, Y stands for the primary variables and X stands for the intervention variables.
  ## Refer to descriptions of simulations related to Table 2 in the main text for details.
  generation <- hub.generation(p,n)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  W <- generation$W
  V <- W%*%solve(diag(p)- U)
  
  out <- topological_order(V)
  Pi <- out$an_mat ## The matrix representing ancestral relationship, 1 if existing, 0 otherwise.
  Phi <- out$in_mat ## The matrix representing invertion relationship, 1 if existing, 0 otherwise.
  Piv <- out$iv_mat ## The matrix representing candidate IV set, 1 if existing, 0 otherwise.
  
  
  ## Compute the evaluation metrics for method without taking consideration of confounders in parameter estimation. 
  ## u1: estimated coefficient matrix of DAG by directly sparse regression.
  u1 <- coef.direct(X,Y,Pi,Phi,tau.list,gamma.list,n.fold)
  res[1:3] <- metrics(U,u1,Pi)
  
  ## Compute the evaluation metrics for GrIVET in parameter estimation. 
  ## u2: estimated coefficient matrix of DAG by GrIVET.
  u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  res[4:6] <- metrics(U,u2,Pi)
  
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 6
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);stat <- simu(i);return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./primary_results/part2/setting1/","stats_hub3.csv"),row.names = FALSE)

## hub graph, n = 500, discrete intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 101;n <- 500;
  tau.list <- seq(0.2,0.3,0.01)
  gamma.list <- seq(0.1,1,0.1)
  n.fold <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(6) # `res` is to store the metrics for the evaluation.
  
  ## Data generation, Y stands for the primary variables and X stands for the intervention variables.
  ## Refer to descriptions of simulations related to Table 2 in the main text for details.
  generation <- hub.generation2(p,n)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  W <- generation$W
  V <- W%*%solve(diag(p)- U)
  
  out <- topological_order(V)
  Pi <- out$an_mat ## The matrix representing ancestral relationship, 1 if existing, 0 otherwise.
  Phi <- out$in_mat ## The matrix representing invertion relationship, 1 if existing, 0 otherwise.
  Piv <- out$iv_mat ## The matrix representing candidate IV set, 1 if existing, 0 otherwise.
  
  
  ## Compute the evaluation metrics for method without taking consideration of confounders in parameter estimation. 
  ## u1: estimated coefficient matrix of DAG by directly sparse regression.
  u1 <- coef.direct(X,Y,Pi,Phi,tau.list,gamma.list,n.fold)
  res[1:3] <- metrics(U,u1,Pi)
  
  ## Compute the evaluation metrics for GrIVET in parameter estimation. 
  ## u2: estimated coefficient matrix of DAG by GrIVET.
  u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  res[4:6] <- metrics(U,u2,Pi)
  
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 6
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);stat <- simu(i);return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./primary_results/part2/setting1/","stats_hub4.csv"),row.names = FALSE)

## hub graph, n = 400, discrete intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 101;n <- 400;
  tau.list <- seq(0.2,0.3,0.01)
  gamma.list <- seq(0.1,1,0.1)
  n.fold <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(6) # `res` is to store the metrics for the evaluation.
  
  ## Data generation, Y stands for the primary variables and X stands for the intervention variables.
  ## Refer to descriptions of simulations related to Table 2 in the main text for details.
  generation <- hub.generation2(p,n)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  W <- generation$W
  V <- W%*%solve(diag(p)- U)
  
  out <- topological_order(V)
  Pi <- out$an_mat ## The matrix representing ancestral relationship, 1 if existing, 0 otherwise.
  Phi <- out$in_mat ## The matrix representing invertion relationship, 1 if existing, 0 otherwise.
  Piv <- out$iv_mat ## The matrix representing candidate IV set, 1 if existing, 0 otherwise.
  
  
  ## Compute the evaluation metrics for method without taking consideration of confounders in parameter estimation. 
  ## u1: estimated coefficient matrix of DAG by directly sparse regression.
  u1 <- coef.direct(X,Y,Pi,Phi,tau.list,gamma.list,n.fold)
  res[1:3] <- metrics(U,u1,Pi)
  
  ## Compute the evaluation metrics for GrIVET in parameter estimation. 
  ## u2: estimated coefficient matrix of DAG by GrIVET.
  u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  res[4:6] <- metrics(U,u2,Pi)
  
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 6
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);stat <- simu(i);return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./primary_results/part2/setting1/","stats_hub5.csv"),row.names = FALSE)


## hub graph, n = 300, discrete intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 101;n <- 300;
  tau.list <- seq(0.2,0.3,0.01)
  gamma.list <- seq(0.1,1,0.1)
  n.fold <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(6) # `res` is to store the metrics for the evaluation.
  
  ## Data generation, Y stands for the primary variables and X stands for the intervention variables.
  ## Refer to descriptions of simulations related to Table 2 in the main text for details.
  generation <- hub.generation2(p,n)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  W <- generation$W
  V <- W%*%solve(diag(p)- U)
  
  out <- topological_order(V)
  Pi <- out$an_mat ## The matrix representing ancestral relationship, 1 if existing, 0 otherwise.
  Phi <- out$in_mat ## The matrix representing invertion relationship, 1 if existing, 0 otherwise.
  Piv <- out$iv_mat ## The matrix representing candidate IV set, 1 if existing, 0 otherwise.
  
  
  ## Compute the evaluation metrics for method without taking consideration of confounders in parameter estimation. 
  ## u1: estimated coefficient matrix of DAG by directly sparse regression.
  u1 <- coef.direct(X,Y,Pi,Phi,tau.list,gamma.list,n.fold)
  res[1:3] <- metrics(U,u1,Pi)
  
  ## Compute the evaluation metrics for GrIVET in parameter estimation. 
  ## u2: estimated coefficient matrix of DAG by GrIVET.
  u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  res[4:6] <- metrics(U,u2,Pi)
  
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 6
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);stat <- simu(i);return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./primary_results/part2/setting1/","stats_hub6.csv"),row.names = FALSE)


stopCluster(cl)