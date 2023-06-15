# This file produces simulations corresponding to Table 1 in the supplementary materials.

## In the following `simu` functions, unless explicitly specified, 
## p: the number of primary variables.
## n: the number of observations.
## tau.list1 & gamma.list1: lists of tunning parameters for TLP regression in matrix V estimation.
## tau.list2 & gamma.list2; lists of tunning parameters for TLP regression in structure learning(penalized coefficient estimation).
## n.fold1 & n.fold2: number of folds for cross-validation in tuning parameters selection, respectively for V estimation and structure learning.
## type: type = 0 controls the simulated graph so that H_0 holds, otherwise H_1 holds.

library(doParallel)
n_cores <- detectCores() - 1
registerDoParallel(cores=n_cores)  
cl <- makeCluster(n_cores, type="FORK")

clusterEvalQ(cl,{source("comparison1_setting2.R",chdir = TRUE)
  library(MASS)})

## random graph, n = 1000, continuous intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 10;n <- 1000;
  tau.list1 <- seq(0.1,0.2,0.01)
  gamma.list1 <- seq(0.05,0.5,0.05)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  n.fold1 <- n.fold2 <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(12) # `res` is to store the metrics for the evaluation of different methods.
  
  ## Data generation, Y stands for the primary variables and X stands for the intervention variables.
  ## Refer to descriptions of simulations related to Table 1 in the supplementary materials for details.
  generation <- random.generation(p,n,U.test,type)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  U <- 1*(U!=0) # record the support of structure for DAG.
  y <- cbind(Y,X)
  
  
  ## Compute the evaluation metrics for RFCI. 
  ## u1: estimated structure of DAG by RFCI.
  out1 <- causal_discovery_RFCI(y)
  u1 <- out1$u
  u1 <- u1[1:p,1:p]
  
  u1 <- PAG2best(U,u1)
  res[1:4] <- metrics(U,u1)
  
  ## Compute the evaluation metrics for LRpS_GES. 
  ## u2: estimated structure of DAG by LRpS_GES.
  out2 <- causal_discovery_LRpS_GES(y)
  u2 <- out2$u
  u2 <- u2[1:p,1:p]
  res[5:8] <- metrics(U,u2)
  
  ## Compute the evaluation metrics for GrIVET. 
  ## u3: estimated structure of DAG by GrIVET.
  u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
  res[9:12] <- metrics(U,u3)
  
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 12
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);stat<-simu(i);return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=12)
write.csv(stats_random,file.path("./primary_results/part1/setting2/","stats_random1.csv"),row.names = FALSE)

## random graph, n = 500, continuous intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 10;n <- 500;
  tau.list1 <- seq(0.1,0.2,0.01)
  gamma.list1 <- seq(0.05,0.5,0.05)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  n.fold1 <- n.fold2 <- 5
  U.test <- matrix(0,p,p) # record the support of structure for DAG.
  type <- 0
  res <- numeric(12) # `res` is to store the metrics for the evaluation of different methods.
  
  ## Data generation, Y stands for the primary variables and X stands for the intervention variables.
  ## Refer to descriptions of simulations related to Table 1 in the supplementary materials for details.
  generation <- random.generation(p,n,U.test,type)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  U <- 1*(U!=0)
  y <- cbind(Y,X)
  
  
  ## Compute the evaluation metrics for RFCI. 
  ## u1: estimated structure of DAG by RFCI.
  out1 <- causal_discovery_RFCI(y)
  u1 <- out1$u
  u1 <- u1[1:p,1:p]
  
  u1 <- PAG2best(U,u1)
  res[1:4] <- metrics(U,u1)
  
  ## Compute the evaluation metrics for LRpS_GES. 
  ## u2: estimated structure of DAG by LRpS_GES.
  out2 <- causal_discovery_LRpS_GES(y)
  u2 <- out2$u
  u2 <- u2[1:p,1:p]
  res[5:8] <- metrics(U,u2)
  
  ## Compute the evaluation metrics for GrIVET. 
  ## u3: estimated structure of DAG by GrIVET.
  u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
  res[9:12] <- metrics(U,u3)
  
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 12
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=12)
write.csv(stats_random,file.path("./primary_results/part1/setting2/","stats_random2.csv"),row.names = FALSE)


## random graph, n = 200, continuous intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 10;n <- 200;
  tau.list1 <- seq(0.1,0.2,0.01)
  gamma.list1 <- seq(0.05,0.5,0.05)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  n.fold1 <- n.fold2 <- 5
  U.test <- matrix(0,p,p) # record the support of structure for DAG.
  type <- 0
  res <- numeric(12) # `res` is to store the metrics for the evaluation of different methods.
  
  ## Data generation, Y stands for the primary variables and X stands for the intervention variables.
  ## Refer to descriptions of simulations related to Table 1 in the supplementary materials for details.
  generation <- random.generation(p,n,U.test,type)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  U <- 1*(U!=0)
  y <- cbind(Y,X)
  
  
  ## Compute the evaluation metrics for RFCI. 
  ## u1: estimated structure of DAG by RFCI.
  out1 <- causal_discovery_RFCI(y)
  u1 <- out1$u
  u1 <- u1[1:p,1:p]
  
  u1 <- PAG2best(U,u1)
  res[1:4] <- metrics(U,u1)
  
  ## Compute the evaluation metrics for LRpS_GES. 
  ## u2: estimated structure of DAG by LRpS_GES.
  out2 <- causal_discovery_LRpS_GES(y)
  u2 <- out2$u
  u2 <- u2[1:p,1:p]
  res[5:8] <- metrics(U,u2)
  
  ## Compute the evaluation metrics for GrIVET. 
  ## u3: estimated structure of DAG by GrIVET.
  u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
  res[9:12] <- metrics(U,u3)
  
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 12
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=12)
write.csv(stats_random,file.path("./primary_results/part1/setting2/","stats_random3.csv"),row.names = FALSE)

## random graph, n = 1000, discrete intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 10;n <- 1000;
  tau.list1 <- seq(0.1,0.2,0.01)
  gamma.list1 <- seq(0.05,0.5,0.05)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  n.fold1 <- n.fold2 <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(12) # `res` is to store the metrics for the evaluation of different methods.
  
  ## Data generation, Y stands for the primary variables and X stands for the intervention variables.
  ## Refer to descriptions of simulations related to Table 1 in the supplementary materials for details.
  generation <- random.generation2(p,n,U.test,type)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  U <- 1*(U!=0) # record the support of structure for DAG.
  y <- cbind(Y,X)
  
  
  ## Compute the evaluation metrics for RFCI. 
  ## u1: estimated structure of DAG by RFCI.
  out1 <- causal_discovery_RFCI(y)
  u1 <- out1$u
  u1 <- u1[1:p,1:p]
  
  u1 <- PAG2best(U,u1)
  res[1:4] <- metrics(U,u1)
  
  ## Compute the evaluation metrics for LRpS_GES. 
  ## u2: estimated structure of DAG by LRpS_GES.
  out2 <- causal_discovery_LRpS_GES(y)
  u2 <- out2$u
  u2 <- u2[1:p,1:p]
  res[5:8] <- metrics(U,u2)
  
  ## Compute the evaluation metrics for GrIVET. 
  ## u3: estimated structure of DAG by GrIVET.
  u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
  res[9:12] <- metrics(U,u3)
  
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 12
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=12)
write.csv(stats_random,file.path("./primary_results/part1/setting2/","stats_random4.csv"),row.names = FALSE)

## random graph, n = 500, discrete intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 10;n <- 500;
  tau.list1 <- seq(0.1,0.2,0.01)
  gamma.list1 <- seq(0.05,0.5,0.05)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  n.fold1 <- n.fold2 <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(12) # `res` is to store the metrics for the evaluation of different methods.
  
  ## Data generation, Y stands for the primary variables and X stands for the intervention variables.
  ## Refer to descriptions of simulations related to Table 1 in the supplementary materials for details.
  generation <- random.generation2(p,n,U.test,type)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  U <- 1*(U!=0) # record the support of structure for DAG.
  y <- cbind(Y,X)
  
  
  ## Compute the evaluation metrics for RFCI. 
  ## u1: estimated structure of DAG by RFCI.
  out1 <- causal_discovery_RFCI(y)
  u1 <- out1$u
  u1 <- u1[1:p,1:p]
  
  u1 <- PAG2best(U,u1)
  res[1:4] <- metrics(U,u1)
  
  ## Compute the evaluation metrics for LRpS_GES. 
  ## u2: estimated structure of DAG by LRpS_GES.
  out2 <- causal_discovery_LRpS_GES(y)
  u2 <- out2$u
  u2 <- u2[1:p,1:p]
  res[5:8] <- metrics(U,u2)
  
  ## Compute the evaluation metrics for GrIVET. 
  ## u3: estimated structure of DAG by GrIVET.
  u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
  res[9:12] <- metrics(U,u3)
  
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 12
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=12)
write.csv(stats_random,file.path("./primary_results/part1/setting2/","stats_random5.csv"),row.names = FALSE)


## random graph, n = 200, discrete intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 10;n <- 200;
  tau.list1 <- seq(0.1,0.2,0.01)
  gamma.list1 <- seq(0.05,0.5,0.05)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  n.fold1 <- n.fold2 <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(12) # `res` is to store the metrics for the evaluation of different methods.
  
  ## Data generation, Y stands for the primary variables and X stands for the intervention variables.
  ## Refer to descriptions of simulations related to Table 1 in the supplementary materials for details.
  generation <- random.generation2(p,n,U.test,type)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  U <- 1*(U!=0) # record the support of structure for DAG.
  y <- cbind(Y,X)
  
  
  ## Compute the evaluation metrics for RFCI. 
  ## u1: estimated structure of DAG by RFCI.
  out1 <- causal_discovery_RFCI(y)
  u1 <- out1$u
  u1 <- u1[1:p,1:p]
  
  u1 <- PAG2best(U,u1)
  res[1:4] <- metrics(U,u1)
  
  ## Compute the evaluation metrics for LRpS_GES. 
  ## u2: estimated structure of DAG by LRpS_GES.
  out2 <- causal_discovery_LRpS_GES(y)
  u2 <- out2$u
  u2 <- u2[1:p,1:p]
  res[5:8] <- metrics(U,u2)
  
  ## Compute the evaluation metrics for GrIVET. 
  ## u3: estimated structure of DAG by GrIVET.
  u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
  res[9:12] <- metrics(U,u3)
  
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 12
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=12)
write.csv(stats_random,file.path("./primary_results/part1/setting2/","stats_random6.csv"),row.names = FALSE)

## hub graph, n = 1000, continuous intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 11;n <- 1000;
  tau.list1 <- seq(0.2,0.3,0.01)
  gamma.list1 <- seq(0.1,1,0.1)
  tau.list2 <- seq(0.2,0.3,0.01)
  gamma.list2 <- seq(0.1,1,0.1)
  n.fold1 <- n.fold2 <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(12) # `res` is to store the metrics for the evaluation of different methods.
  
  ## Data generation, Y stands for the primary variables and X stands for the intervention variables.
  ## Refer to descriptions of simulations related to Table 1 in the supplementary materials for details.
  generation <- hub.generation(p,n)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  U <- 1*(U!=0) # record the support of structure for DAG.
  y <- cbind(Y,X)
  
  
  ## Compute the evaluation metrics for RFCI. 
  ## u1: estimated structure of DAG by RFCI.
  out1 <- causal_discovery_RFCI(y)
  u1 <- out1$u
  u1 <- u1[1:p,1:p]
  
  u1 <- PAG2best(U,u1)
  res[1:4] <- metrics(U,u1)
  
  ## Compute the evaluation metrics for LRpS_GES. 
  ## u2: estimated structure of DAG by LRpS_GES.
  out2 <- causal_discovery_LRpS_GES(y)
  u2 <- out2$u
  u2 <- u2[1:p,1:p]
  res[5:8] <- metrics(U,u2)
  
  ## Compute the evaluation metrics for GrIVET. 
  ## u3: estimated structure of DAG by GrIVET.
  u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
  res[9:12] <- metrics(U,u3)
  
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 12
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=12)
write.csv(stats_random,file.path("./primary_results/part1/setting2/","stats_hub1.csv"),row.names = FALSE)

## hub graph, n = 500, continuous intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 11;n <- 500;
  tau.list1 <- seq(0.2,0.3,0.01)
  gamma.list1 <- seq(0.1,1,0.1)
  tau.list2 <- seq(0.2,0.3,0.01)
  gamma.list2 <- seq(0.1,1,0.1)
  n.fold1 <- n.fold2 <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(12) # `res` is to store the metrics for the evaluation of different methods.
  
  ## Data generation, Y stands for the primary variables and X stands for the intervention variables.
  ## Refer to descriptions of simulations related to Table 1 in the supplementary materials for details.
  generation <- hub.generation(p,n)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  U <- 1*(U!=0) # record the support of structure for DAG.
  y <- cbind(Y,X)
  
  
  ## Compute the evaluation metrics for RFCI. 
  ## u1: estimated structure of DAG by RFCI.
  out1 <- causal_discovery_RFCI(y)
  u1 <- out1$u
  u1 <- u1[1:p,1:p]
  
  u1 <- PAG2best(U,u1)
  res[1:4] <- metrics(U,u1)
  
  ## Compute the evaluation metrics for LRpS_GES. 
  ## u2: estimated structure of DAG by LRpS_GES.
  out2 <- causal_discovery_LRpS_GES(y)
  u2 <- out2$u
  u2 <- u2[1:p,1:p]
  res[5:8] <- metrics(U,u2)
  
  ## Compute the evaluation metrics for GrIVET. 
  ## u3: estimated structure of DAG by GrIVET.
  u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
  res[9:12] <- metrics(U,u3)
  
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 12
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=12)
write.csv(stats_random,file.path("./primary_results/part1/setting2/","stats_hub2.csv"),row.names = FALSE)


## hub graph, n = 200, continuous intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 11;n <- 200;
  tau.list1 <- seq(0.2,0.3,0.01)
  gamma.list1 <- seq(0.1,1,0.1)
  tau.list2 <- seq(0.2,0.3,0.01)
  gamma.list2 <- seq(0.1,1,0.1)
  n.fold1 <- n.fold2 <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(12) # `res` is to store the metrics for the evaluation of different methods.
  
  ## Data generation, Y stands for the primary variables and X stands for the intervention variables.
  ## Refer to descriptions of simulations related to Table 1 in the supplementary materials for details.
  generation <- hub.generation(p,n)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  U <- 1*(U!=0) # record the support of structure for DAG.
  y <- cbind(Y,X)
  
  
  ## Compute the evaluation metrics for RFCI. 
  ## u1: estimated structure of DAG by RFCI.
  out1 <- causal_discovery_RFCI(y)
  u1 <- out1$u
  u1 <- u1[1:p,1:p]
  
  u1 <- PAG2best(U,u1)
  res[1:4] <- metrics(U,u1)
  
  ## Compute the evaluation metrics for LRpS_GES. 
  ## u2: estimated structure of DAG by LRpS_GES.
  out2 <- causal_discovery_LRpS_GES(y)
  u2 <- out2$u
  u2 <- u2[1:p,1:p]
  res[5:8] <- metrics(U,u2)
  
  ## Compute the evaluation metrics for GrIVET. 
  ## u3: estimated structure of DAG by GrIVET.
  u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
  res[9:12] <- metrics(U,u3)
  
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 12
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=12)
write.csv(stats_random,file.path("./primary_results/part1/setting2/","stats_hub3.csv"),row.names = FALSE)

## hub graph, n = 1000, discrete intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 11;n <- 1000;
  tau.list1 <- seq(0.2,0.3,0.01)
  gamma.list1 <- seq(0.1,1,0.1)
  tau.list2 <- seq(0.2,0.3,0.01)
  gamma.list2 <- seq(0.1,1,0.1)
  n.fold1 <- n.fold2 <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(12) # `res` is to store the metrics for the evaluation of different methods.
  
  ## Data generation, Y stands for the primary variables and X stands for the intervention variables.
  ## Refer to descriptions of simulations related to Table 1 in the supplementary materials for details.
  generation <- hub.generation2(p,n)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  U <- 1*(U!=0) # record the support of structure for DAG.
  y <- cbind(Y,X)
  
  
  ## Compute the evaluation metrics for RFCI. 
  ## u1: estimated structure of DAG by RFCI.
  out1 <- causal_discovery_RFCI(y)
  u1 <- out1$u
  u1 <- u1[1:p,1:p]
  
  u1 <- PAG2best(U,u1)
  res[1:4] <- metrics(U,u1)
  
  ## Compute the evaluation metrics for LRpS_GES. 
  ## u2: estimated structure of DAG by LRpS_GES.
  out2 <- causal_discovery_LRpS_GES(y)
  u2 <- out2$u
  u2 <- u2[1:p,1:p]
  res[5:8] <- metrics(U,u2)
  
  ## Compute the evaluation metrics for GrIVET. 
  ## u3: estimated structure of DAG by GrIVET.
  u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
  res[9:12] <- metrics(U,u3)
  
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 12
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=12)
write.csv(stats_random,file.path("./primary_results/part1/setting2/","stats_hub4.csv"),row.names = FALSE)

## hub graph, n = 500, discrete intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 11;n <- 500;
  tau.list1 <- seq(0.2,0.3,0.01)
  gamma.list1 <- seq(0.1,1,0.1)
  tau.list2 <- seq(0.2,0.3,0.01)
  gamma.list2 <- seq(0.1,1,0.1)
  n.fold1 <- n.fold2 <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(12) # `res` is to store the metrics for the evaluation of different methods.
  
  ## Data generation, Y stands for the primary variables and X stands for the intervention variables.
  ## Refer to descriptions of simulations related to Table 1 in the supplementary materials for details.
  generation <- hub.generation2(p,n)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  U <- 1*(U!=0) # record the support of structure for DAG.
  y <- cbind(Y,X)
  
  
  ## Compute the evaluation metrics for RFCI. 
  ## u1: estimated structure of DAG by RFCI.
  out1 <- causal_discovery_RFCI(y)
  u1 <- out1$u
  u1 <- u1[1:p,1:p]
  
  u1 <- PAG2best(U,u1)
  res[1:4] <- metrics(U,u1)
  
  ## Compute the evaluation metrics for LRpS_GES. 
  ## u2: estimated structure of DAG by LRpS_GES.
  out2 <- causal_discovery_LRpS_GES(y)
  u2 <- out2$u
  u2 <- u2[1:p,1:p]
  res[5:8] <- metrics(U,u2)
  
  ## Compute the evaluation metrics for GrIVET. 
  ## u3: estimated structure of DAG by GrIVET.
  u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
  res[9:12] <- metrics(U,u3)
  
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 12
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=12)
write.csv(stats_random,file.path("./primary_results/part1/setting2/","stats_hub5.csv"),row.names = FALSE)


## hub graph, n = 200, discrete intervention variables
simu <- function(i){
  ## Set the parameters for the simulation, refer to the start of the file for the details.
  p <- 11;n <- 200;
  tau.list1 <- seq(0.2,0.3,0.01)
  gamma.list1 <- seq(0.1,1,0.1)
  tau.list2 <- seq(0.2,0.3,0.01)
  gamma.list2 <- seq(0.1,1,0.1)
  n.fold1 <- n.fold2 <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(12) # `res` is to store the metrics for the evaluation of different methods.
  
  ## Data generation, Y stands for the primary variables and X stands for the intervention variables.
  ## Refer to descriptions of simulations related to Table 1 in the supplementary materials for details.
  generation <- hub.generation2(p,n)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  U <- 1*(U!=0) # record the support of structure for DAG.
  y <- cbind(Y,X)
  
  
  ## Compute the evaluation metrics for RFCI. 
  ## u1: estimated structure of DAG by RFCI.
  out1 <- causal_discovery_RFCI(y)
  u1 <- out1$u
  u1 <- u1[1:p,1:p]
  
  u1 <- PAG2best(U,u1)
  res[1:4] <- metrics(U,u1)
  
  ## Compute the evaluation metrics for LRpS_GES. 
  ## u2: estimated structure of DAG by LRpS_GES.
  out2 <- causal_discovery_LRpS_GES(y)
  u2 <- out2$u
  u2 <- u2[1:p,1:p]
  res[5:8] <- metrics(U,u2)
  
  ## Compute the evaluation metrics for GrIVET. 
  ## u3: estimated structure of DAG by GrIVET.
  u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
  res[9:12] <- metrics(U,u3)
  
  return(res)
}

## repeat the simulations for 1000 times and store the results of evaluations.
len.stat <- 12
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=12)
write.csv(stats_random,file.path("./primary_results/part1/setting2/","stats_hub6.csv"),row.names = FALSE)

stopCluster(cl)