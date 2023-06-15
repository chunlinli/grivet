## This file conduct simulations for the hub graph part in table 3 of the supplementary materials.
library(doParallel)
n_cores <- detectCores() - 1
registerDoParallel(cores=n_cores)  
cl <- makeCluster(n_cores, type="FORK")

clusterEvalQ(cl,{library(grivet)
  library(mvtnorm)
  library(MASS)
  library(clusterGeneration)
  library(mnormt)})

## This function computes the test statistic proposed in the GrIVET paper.

## Arguments:
## X: n*q data matrix for intervention variables.
## Y: n*p data matrix for primary variables.
## tau.list1 & gamma.list1: lists of tuning parameters for TLP regression in matrix V estimation.
## tau.list2 & gamma.list2; lists of tuning parameters for TLP regression in parameter estimation.
## tau.list3 & gamma.list3; lists of tuning parameters for TLP regression in support recovery of precision matrix.
## n.fold1 & n.fold2 & n.fold3: number of folds for cross-validation in tunning parameters selection, respectively for V estimation, parameter estimation, precision matrix support recovery.
## U.test.list: a list of matrices representing edges to be tested.
## Sigma_true: the true covariance matrix in DAG.
## max.it & tol: the number of iterations and tolerance in precision matrix refiting.

## Returns:
## statistic: the proposed statistic

simu.single.2lr <- function(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list,Sigma_true,max.it=100000,tol=1e-7){
  ## Split sample
  n <- nrow(X)
  X2 <- X[(n-199):n,]
  Y2 <- Y[(n-199):n,]
  X <- X[1:(n-200),]
  Y <- Y[1:(n-200),]
  
  # stage 1: discovery
  result.stage1_1 <- cv.intdag.pmle.diff.aic(X,Y,tau.list1,gamma.list1,n.fold1)
  V.es <- result.stage1_1$V
  result.stage1_2 <- topological_order(V.es)
  Pi_es <- result.stage1_2$an_mat
  Phi_es <- result.stage1_2$in_mat
  Piv_es <- result.stage1_2$iv_mat
  
  ## Obtain covariance matrix estimation with respect to Pi_es
  result.stage2 <- cv.intdag.coe(X,Y,Pi_es,Phi_es,Piv_es,tau.list2,gamma.list2,n.fold2)
  Z <- Y - Y%*%result.stage2$U - X%*%result.stage2$W
  MB.support <- cv.MB_Union(Z,tau.list3,gamma.list3,n.fold3)$S
  Sigma_es <- result.stage2$Sigma
  S_es <- precision_refit(Sigma_es,MB.support,max.it,tol)
  
  
  ## Compute 2lr using the whole graph
  statistic <- numeric(length(U.test.list))
  for(k in 1:length(U.test.list)){
    U.test <- U.test.list[[k]]
    ## combine U.test with Pi_es
    U1 <- 1*((Pi_es+U.test)>0)
    if (is.acyclic(U1)==0){
      statistic = 0
      return(statistic)
    }else{
      U0 <- U1-U.test
    }
    try(result.stage3 <- intdag.2lr(X2,Y2,U0,Phi_es,U1,Phi_es,S_es)$statistic,result.stage3 <- -1)
    statistic[k] <- result.stage3
  }
  
  return(statistic)
}

## This function judeges if the given `U` encodes a directed acyclic graph.

## Arguments:
## U: the causal effect matrix.

## Returns:
## flag: the encoded graph has a cycle if flag = 0, acyclic if flag = 1.
is.acyclic <- function(U){
  flag <- 1
  while (sum(U)>0){
    if (min(colSums(U))>0){
      flag <- 0
      break
    }
    idx <- which(colSums(U)!=0)
    U <- U[idx,idx,drop=FALSE]
  }
  return(flag)
}

## This function generate a hub graph with continuous intervention variables.

## Arguments: 
## p: the number of primary variables.
## n: the number of observations.
## U.test: a p*p matrix, 1 representing tested edges, 0 otherwise.
## type: type = 0 controls the simulated graph so that H_0 holds, otherwise H_1 holds.

## Returns:
## X: a n*q matrix for intervention variables.
## Y: a n*p matrix for primary variables.
## U: a p*P matrix representing the causal effect matrix.
## W: a q*p matrix representing the interventional effect matrix.
## Sigma: a p*p matrix representing the covariance matrix of residuals.
hub.generation <- function(p,n){
  if (p%%2==0)
    stop("wrong p is given")
  q <- 2*p+(p-1)/2
  U <- matrix(0,p,p);U[1,2:p] <- 1;
  W1 <- W2 <- diag(1,p,p); W3 <- matrix(0,q-2*p,p); 
  for (i in 1:(q-2*p)){
    W3[i,2*i]<-1
    W3[i,(2*i+1)]<-1
  }
  W <- rbind(W1,W2,W3)
  U <- U*matrix(runif(p*p,min=0.8,max=1.2)*sample(c(-1,1),p*p,replace = TRUE),p,p)
  W <- W*matrix(runif(q*p,min=0.8,max=1.2)*sample(c(-1,1),q*p,replace = TRUE),q,p)
  vari <- diag(runif(p,min=0.4,max=0.6))
  sigma <- vari%*%diag(p)%*%vari
  e <- mvrnorm(n,rep(0,p),sigma)
  r <- (p-1)/2
  V <- matrix(0,r,p)
  for (i in 1:r){
    V[i,2*i] <- 1
    V[i,2*i+1] <- 1
  }
  V[1,1] <- 1
  V <- V*matrix(runif(r*p,min=0.4,max=0.6)*sample(c(-1,1),r*p,replace = TRUE),r,p)
  Z <- matrix(rnorm(r*n),ncol=r)
  X <- matrix(rnorm(q*n),ncol=q); Y <- matrix(0,n,p); Y <- (X%*%W+e+Z%*%V)%*%solve(diag(1,p)-U);
  Sigma <- sigma+t(V)%*%V
  return(list(X=X,Y=Y,U=U,W=W,Sigma=Sigma))
}

## This function generate a hub graph with discrete intervention variables.

## Arguments: 
## p: the number of primary variables.
## n: the number of observations.
## U.test: a p*p matrix, 1 representing tested edges, 0 otherwise.
## type: type = 0 controls the simulated graph so that H_0 holds, otherwise H_1 holds.

## Returns:
## X: a n*q matrix for intervention variables.
## Y: a n*p matrix for primary variables.
## U: a p*P matrix representing the causal effect matrix.
## W: a q*p matrix representing the interventional effect matrix.
## Sigma: a p*p matrix representing the covariance matrix of residuals.
hub.generation2 <- function(p,n){
  if (p%%2==0)
    stop("wrong p is given")
  q <- 2*p+(p-1)/2
  U <- matrix(0,p,p);U[1,2:p] <- 1;
  W1 <- W2 <- diag(1,p,p); W3 <- matrix(0,q-2*p,p); 
  for (i in 1:(q-2*p)){
    W3[i,2*i]<-1
    W3[i,(2*i+1)]<-1
  }
  W <- rbind(W1,W2,W3)
  U <- U*matrix(runif(p*p,min=0.8,max=1.2)*sample(c(-1,1),p*p,replace = TRUE),p,p)
  W <- W*matrix(runif(q*p,min=0.8,max=1.2)*sample(c(-1,1),q*p,replace = TRUE),q,p)
  vari <- diag(runif(p,min=0.4,max=0.6))
  sigma <- vari%*%diag(p)%*%vari
  e <- mvrnorm(n,rep(0,p),sigma)
  r <- (p-1)/2
  V <- matrix(0,r,p)
  for (i in 1:r){
    V[i,2*i] <- 1
    V[i,2*i+1] <- 1
  }
  V[1,1] <- 1
  V <- V*matrix(runif(r*p,min=0.4,max=0.6)*sample(c(-1,1),r*p,replace = TRUE),r,p)
  Z <- matrix(rnorm(r*n),ncol=r)
  X <- matrix(sample(c(-1,1),q*n,replace = TRUE),ncol=q); Y <- matrix(0,n,p); Y <- (X%*%W+e+Z%*%V)%*%solve(diag(1,p)-U);
  Sigma <- sigma+t(V)%*%V
  return(list(X=X,Y=Y,U=U,W=W,Sigma=Sigma))
}

## This function generates hub graphs with continuous interventions and compute statistics with respect to the tested edges.

## Arguments: 
## times: the repeated times.
## p: the number of primary variables.
## n: the number of observations.
## tau.list1 & gamma.list1: lists of tuning parameters for TLP regression in matrix V estimation.
## tau.list2 & gamma.list2; lists of tuning parameters for TLP regression in parameter estimation.
## tau.list3 & gamma.list3; lists of tuning parameters for TLP regression in support recovery of precision matrix.
## n.fold1 & n.fold2 & n.fold3: number of folds for cross-validation in tunning parameters selection, respectively for V estimation, parameter estimation, precision matrix support recovery.
## U.test.list: a list of matrices representing edges to be tested.
## type: type = 0 controls the simulated graph so that H_0 holds, otherwise H_1 holds.

## Returns:
## stats: the computed test statistics.

simu.hub <- function(times,p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list){
  stats <- matrix(0,times,length(U.test.list))
  for (i in 1:times){
    generation <- hub.generation(p,n)
    X <- generation$X
    Y <- generation$Y
    Sigma_true <- generation$Sigma
    stats[i,] <- simu.single.2lr(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list,Sigma_true)
  }
  return(stats)
}

## This function generates hub graphs with discrete interventions and compute statistics with respect to the tested edges.

## Arguments: 
## times: the repeated times.
## p: the number of primary variables.
## n: the number of observations.
## tau.list1 & gamma.list1: lists of tuning parameters for TLP regression in matrix V estimation.
## tau.list2 & gamma.list2; lists of tuning parameters for TLP regression in parameter estimation.
## tau.list3 & gamma.list3; lists of tuning parameters for TLP regression in support recovery of precision matrix.
## n.fold1 & n.fold2 & n.fold3: number of folds for cross-validation in tunning parameters selection, respectively for V estimation, parameter estimation, precision matrix support recovery.
## U.test.list: a list of matrices representing edges to be tested.
## type: type = 0 controls the simulated graph so that H_0 holds, otherwise H_1 holds.

## Returns:
## stats: the computed test statistics.
simu.hub2 <- function(times,p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list){
  stats <- matrix(0,times,length(U.test.list))
  for (i in 1:times){
    generation <- hub.generation2(p,n)
    X <- generation$X
    Y <- generation$Y
    Sigma_true <- generation$Sigma
    stats[i,] <- simu.single.2lr(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list,Sigma_true)
  }
  return(stats)
}

clusterExport(cl,c("simu.single.2lr","is.acyclic","hub.generation","hub.generation2","simu.hub","simu.hub2"))

## hub graph, n = 1000, m = 200, continuous intervention variables
simu <- function(i){
  ## set the parameters, see arguments in function "simu.hub" for details.
  times <- 1;p <- 11;n <- 1200;
  tau.list1 <- seq(0.2,0.3,0.01)
  gamma.list1 <- seq(0.1,1,0.1)
  tau.list2 <- seq(0.2,0.3,0.01)
  gamma.list2 <- seq(0.1,1,0.1)
  tau.list3 <- seq(0.05,0.1,0.01)
  gamma.list3 <- seq(0.005,0.1,0.005)
  n.fold1 <- n.fold2 <- n.fold3 <- 5
  
  ## create two lists of matrices representing the edges to be tested for different testing types.
  U.test.list <- vector("list",6)
  U.test1 <- matrix(0,p,p)
  U.test1[2,3] <- 1;
  U.test2 <- matrix(0,p,p)
  U.test2[2,3] <- 1; U.test2[3,4] <- 1; U.test2[4,5] <- 1;
  U.test3 <- matrix(0,p,p)
  U.test3[2,3] <- 1; U.test3[3,4] <- 1; U.test3[4,5] <- 1; U.test3[5,6] <- 1; U.test3[6,7] <- 1;
  U.test4 <- matrix(0,p,p)
  U.test4[1,2] <- 1;
  U.test5 <- matrix(0,p,p)
  U.test5[1,2] <- 1; U.test5[1,3] <- 1; U.test5[1,4] <- 1;
  U.test6 <- matrix(0,p,p)
  U.test6[1,2] <- 1; U.test6[1,3] <- 1; U.test6[1,4] <- 1; U.test6[1,5] <- 1; U.test6[1,6] <- 1;
  U.test.list[[1]] <- U.test1; U.test.list[[2]] <- U.test2; U.test.list[[3]] <- U.test3;
  U.test.list[[4]] <- U.test4; U.test.list[[5]] <- U.test5; U.test.list[[6]] <- U.test6;
  
  ## compute the test statistics with respect to different numbers of edges to be tested and different testing type.
  stats.hub <- simu.hub(times,p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list)
  return(stats.hub)
}

## repeat the simulations for 1000 times and store the results of test statistics.
len.test <- 6
clusterExport(cl,c("simu","len.test"))
stats.hub <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.test));return(stat)})
stats_hub <- matrix(unlist(stats.hub),byrow = TRUE,ncol=6)
write.csv(stats_hub,file.path("./primary_results/part3/setting2/","stats_hub1.csv"),row.names = FALSE)

## hub graph, n = 500, m = 200, continuous intervention variables
simu <- function(i){
  ## set the parameters, see arguments in function "simu.hub" for details.
  times <- 1;p <- 11;n <- 700;
  tau.list1 <- seq(0.2,0.3,0.01)
  gamma.list1 <- seq(0.1,1,0.1)
  tau.list2 <- seq(0.2,0.3,0.01)
  gamma.list2 <- seq(0.1,1,0.1)
  tau.list3 <- seq(0.05,0.1,0.01)
  gamma.list3 <- seq(0.005,0.1,0.005)
  n.fold1 <- n.fold2 <- n.fold3 <- 5
  
  ## create two lists of matrices representing the edges to be tested for different testing types.
  U.test.list <- vector("list",6)
  U.test1 <- matrix(0,p,p)
  U.test1[2,3] <- 1;
  U.test2 <- matrix(0,p,p)
  U.test2[2,3] <- 1; U.test2[3,4] <- 1; U.test2[4,5] <- 1;
  U.test3 <- matrix(0,p,p)
  U.test3[2,3] <- 1; U.test3[3,4] <- 1; U.test3[4,5] <- 1; U.test3[5,6] <- 1; U.test3[6,7] <- 1;
  U.test4 <- matrix(0,p,p)
  U.test4[1,2] <- 1;
  U.test5 <- matrix(0,p,p)
  U.test5[1,2] <- 1; U.test5[1,3] <- 1; U.test5[1,4] <- 1;
  U.test6 <- matrix(0,p,p)
  U.test6[1,2] <- 1; U.test6[1,3] <- 1; U.test6[1,4] <- 1; U.test6[1,5] <- 1; U.test6[1,6] <- 1;
  U.test.list[[1]] <- U.test1; U.test.list[[2]] <- U.test2; U.test.list[[3]] <- U.test3;
  U.test.list[[4]] <- U.test4; U.test.list[[5]] <- U.test5; U.test.list[[6]] <- U.test6;
  
  ## compute the test statistics with respect to different numbers of edges to be tested and different testing type.
  stats.hub <- simu.hub(times,p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list)
  return(stats.hub)
}

## repeat the simulations for 1000 times and store the results of test statistics.
len.test <- 6
clusterExport(cl,c("simu","len.test"))
stats.hub <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.test));return(stat)})
stats_hub <- matrix(unlist(stats.hub),byrow = TRUE,ncol=6)
write.csv(stats_hub,file.path("./primary_results/part3/setting2/","stats_hub2.csv"),row.names = FALSE)


## hub graph, n = 200, m = 200, continuous intervention variables
simu <- function(i){
  ## set the parameters, see arguments in function "simu.hub" for details.
  times <- 1;p <- 11;n <- 400;
  tau.list1 <- seq(0.2,0.3,0.01)
  gamma.list1 <- seq(0.1,1,0.1)
  tau.list2 <- seq(0.2,0.3,0.01)
  gamma.list2 <- seq(0.1,1,0.1)
  tau.list3 <- seq(0.05,0.1,0.01)
  gamma.list3 <- seq(0.005,0.1,0.005)
  n.fold1 <- n.fold2 <- n.fold3 <- 5
  
  ## create two lists of matrices representing the edges to be tested for different testing types.
  U.test.list <- vector("list",6)
  U.test1 <- matrix(0,p,p)
  U.test1[2,3] <- 1;
  U.test2 <- matrix(0,p,p)
  U.test2[2,3] <- 1; U.test2[3,4] <- 1; U.test2[4,5] <- 1;
  U.test3 <- matrix(0,p,p)
  U.test3[2,3] <- 1; U.test3[3,4] <- 1; U.test3[4,5] <- 1; U.test3[5,6] <- 1; U.test3[6,7] <- 1;
  U.test4 <- matrix(0,p,p)
  U.test4[1,2] <- 1;
  U.test5 <- matrix(0,p,p)
  U.test5[1,2] <- 1; U.test5[1,3] <- 1; U.test5[1,4] <- 1;
  U.test6 <- matrix(0,p,p)
  U.test6[1,2] <- 1; U.test6[1,3] <- 1; U.test6[1,4] <- 1; U.test6[1,5] <- 1; U.test6[1,6] <- 1;
  U.test.list[[1]] <- U.test1; U.test.list[[2]] <- U.test2; U.test.list[[3]] <- U.test3;
  U.test.list[[4]] <- U.test4; U.test.list[[5]] <- U.test5; U.test.list[[6]] <- U.test6;
  
  ## compute the test statistics with respect to different numbers of edges to be tested and different testing type.
  stats.hub <- simu.hub(times,p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list)
  return(stats.hub)
}

## repeat the simulations for 1000 times and store the results of test statistics.
len.test <- 6
clusterExport(cl,c("simu","len.test"))
stats.hub <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.test));return(stat)})
stats_hub <- matrix(unlist(stats.hub),byrow = TRUE,ncol=6)
write.csv(stats_hub,file.path("./primary_results/part3/setting2/","stats_hub3.csv"),row.names = FALSE)

## hub graph, n = 1000, m = 200, discrete intervention variables
simu <- function(i){
  ## set the parameters, see arguments in function "simu.hub2" for details.
  times <- 1;p <- 11;n <- 1200;
  tau.list1 <- seq(0.2,0.3,0.01)
  gamma.list1 <- seq(0.1,1,0.1)
  tau.list2 <- seq(0.2,0.3,0.01)
  gamma.list2 <- seq(0.1,1,0.1)
  tau.list3 <- seq(0.05,0.1,0.01)
  gamma.list3 <- seq(0.005,0.1,0.005)
  n.fold1 <- n.fold2 <- n.fold3 <- 5
  
  ## create two lists of matrices representing the edges to be tested for different testing types.
  U.test.list <- vector("list",6)
  U.test1 <- matrix(0,p,p)
  U.test1[2,3] <- 1;
  U.test2 <- matrix(0,p,p)
  U.test2[2,3] <- 1; U.test2[3,4] <- 1; U.test2[4,5] <- 1;
  U.test3 <- matrix(0,p,p)
  U.test3[2,3] <- 1; U.test3[3,4] <- 1; U.test3[4,5] <- 1; U.test3[5,6] <- 1; U.test3[6,7] <- 1;
  U.test4 <- matrix(0,p,p)
  U.test4[1,2] <- 1;
  U.test5 <- matrix(0,p,p)
  U.test5[1,2] <- 1; U.test5[1,3] <- 1; U.test5[1,4] <- 1;
  U.test6 <- matrix(0,p,p)
  U.test6[1,2] <- 1; U.test6[1,3] <- 1; U.test6[1,4] <- 1; U.test6[1,5] <- 1; U.test6[1,6] <- 1;
  U.test.list[[1]] <- U.test1; U.test.list[[2]] <- U.test2; U.test.list[[3]] <- U.test3;
  U.test.list[[4]] <- U.test4; U.test.list[[5]] <- U.test5; U.test.list[[6]] <- U.test6;
  
  ## compute the test statistics with respect to different numbers of edges to be tested and different testing type.
  stats.hub <- simu.hub2(times,p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list)
  return(stats.hub)
}

## repeat the simulations for 1000 times and store the results of test statistics.
len.test <- 6
clusterExport(cl,c("simu","len.test"))
stats.hub <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.test));return(stat)})
stats_hub <- matrix(unlist(stats.hub),byrow = TRUE,ncol=6)
write.csv(stats_hub,file.path("./primary_results/part3/setting2/","stats_hub4.csv"),row.names = FALSE)

## hub graph, n = 500, m = 200, discrete intervention variables
simu <- function(i){
  ## set the parameters, see arguments in function "simu.hub2" for details.
  times <- 1;p <- 11;n <- 700;
  tau.list1 <- seq(0.2,0.3,0.01)
  gamma.list1 <- seq(0.1,1,0.1)
  tau.list2 <- seq(0.2,0.3,0.01)
  gamma.list2 <- seq(0.1,1,0.1)
  tau.list3 <- seq(0.05,0.1,0.01)
  gamma.list3 <- seq(0.005,0.1,0.005)
  n.fold1 <- n.fold2 <- n.fold3 <- 5
  
  ## create two lists of matrices representing the edges to be tested for different testing types.
  U.test.list <- vector("list",6)
  U.test1 <- matrix(0,p,p)
  U.test1[2,3] <- 1;
  U.test2 <- matrix(0,p,p)
  U.test2[2,3] <- 1; U.test2[3,4] <- 1; U.test2[4,5] <- 1;
  U.test3 <- matrix(0,p,p)
  U.test3[2,3] <- 1; U.test3[3,4] <- 1; U.test3[4,5] <- 1; U.test3[5,6] <- 1; U.test3[6,7] <- 1;
  U.test4 <- matrix(0,p,p)
  U.test4[1,2] <- 1;
  U.test5 <- matrix(0,p,p)
  U.test5[1,2] <- 1; U.test5[1,3] <- 1; U.test5[1,4] <- 1;
  U.test6 <- matrix(0,p,p)
  U.test6[1,2] <- 1; U.test6[1,3] <- 1; U.test6[1,4] <- 1; U.test6[1,5] <- 1; U.test6[1,6] <- 1;
  U.test.list[[1]] <- U.test1; U.test.list[[2]] <- U.test2; U.test.list[[3]] <- U.test3;
  U.test.list[[4]] <- U.test4; U.test.list[[5]] <- U.test5; U.test.list[[6]] <- U.test6;
  
  ## compute the test statistics with respect to different numbers of edges to be tested and different testing type.
  stats.hub <- simu.hub2(times,p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list)
  return(stats.hub)
}

## repeat the simulations for 1000 times and store the results of test statistics.
len.test <- 6
clusterExport(cl,c("simu","len.test"))
stats.hub <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.test));return(stat)})
stats_hub <- matrix(unlist(stats.hub),byrow = TRUE,ncol=6)
write.csv(stats_hub,file.path("./primary_results/part3/setting2/","stats_hub5.csv"),row.names = FALSE)


## hub graph, n = 200, m = 200, discrete intervention variables
simu <- function(i){
  ## set the parameters, see arguments in function "simu.hub2" for details.
  times <- 1;p <- 11;n <- 400;
  tau.list1 <- seq(0.2,0.3,0.01)
  gamma.list1 <- seq(0.1,1,0.1)
  tau.list2 <- seq(0.2,0.3,0.01)
  gamma.list2 <- seq(0.1,1,0.1)
  tau.list3 <- seq(0.05,0.1,0.01)
  gamma.list3 <- seq(0.005,0.1,0.005)
  n.fold1 <- n.fold2 <- n.fold3 <- 5
  
  ## create two lists of matrices representing the edges to be tested for different testing types.
  U.test.list <- vector("list",6)
  U.test1 <- matrix(0,p,p)
  U.test1[2,3] <- 1;
  U.test2 <- matrix(0,p,p)
  U.test2[2,3] <- 1; U.test2[3,4] <- 1; U.test2[4,5] <- 1;
  U.test3 <- matrix(0,p,p)
  U.test3[2,3] <- 1; U.test3[3,4] <- 1; U.test3[4,5] <- 1; U.test3[5,6] <- 1; U.test3[6,7] <- 1;
  U.test4 <- matrix(0,p,p)
  U.test4[1,2] <- 1;
  U.test5 <- matrix(0,p,p)
  U.test5[1,2] <- 1; U.test5[1,3] <- 1; U.test5[1,4] <- 1;
  U.test6 <- matrix(0,p,p)
  U.test6[1,2] <- 1; U.test6[1,3] <- 1; U.test6[1,4] <- 1; U.test6[1,5] <- 1; U.test6[1,6] <- 1;
  U.test.list[[1]] <- U.test1; U.test.list[[2]] <- U.test2; U.test.list[[3]] <- U.test3;
  U.test.list[[4]] <- U.test4; U.test.list[[5]] <- U.test5; U.test.list[[6]] <- U.test6;
  
  ## compute the test statistics with respect to different numbers of edges to be tested and different testing type.
  stats.hub <- simu.hub2(times,p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list)
  return(stats.hub)
}

## repeat the simulations for 1000 times and store the results of test statistics.
len.test <- 6
clusterExport(cl,c("simu","len.test"))
stats.hub <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.test));return(stat)})
stats_hub <- matrix(unlist(stats.hub),byrow = TRUE,ncol=6)
write.csv(stats_hub,file.path("./primary_results/part3/setting2/","stats_hub6.csv"),row.names = FALSE)

stopCluster(cl)


