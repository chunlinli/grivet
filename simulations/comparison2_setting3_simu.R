library(doParallel)
n_cores <- detectCores() - 1
registerDoParallel(cores=n_cores)  
cl <- makeCluster(n_cores, type="FORK")

clusterEvalQ(cl,{source("comparison2_setting3.R",chdir = TRUE)
  library(mvtnorm)
  library(MASS)})

## random graph 1
simu <- function(i){
  p <- 10;n <- 1000;
  tau.list <- seq(0.1,0.2,0.01)
  gamma.list <- seq(0.05,0.5,0.05)
  n.fold <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(6)
  
  ## Data generation
  generation <- random.generation(p,n,U.test,type)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  W <- generation$W
  V <- W%*%solve(diag(p)- U)
  
  out <- topological_order(V)
  Pi <- out$an_mat
  Phi <- out$in_mat
  Piv <- out$iv_mat
  
  
  ## For direct computing
  u1 <- coef.direct(X,Y,Pi,Phi,tau.list,gamma.list,n.fold)
  res[1:3] <- metrics(U,u1,Pi)
  
  ## For proposed
  u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  res[4:6] <- metrics(U,u2,Pi)
  
  return(res)
}

len.stat <- 6
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./primary_results/part2/setting3/","stats_random1.csv"),row.names = FALSE)


## random graph 2
simu <- function(i){
  p <- 10;n <- 500;
  tau.list <- seq(0.1,0.2,0.01)
  gamma.list <- seq(0.05,0.5,0.05)
  n.fold <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(6)
  
  ## Data generation
  generation <- random.generation(p,n,U.test,type)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  W <- generation$W
  V <- W%*%solve(diag(p)- U)
  
  out <- topological_order(V)
  Pi <- out$an_mat
  Phi <- out$in_mat
  Piv <- out$iv_mat
  
  
  ## For direct computing
  u1 <- coef.direct(X,Y,Pi,Phi,tau.list,gamma.list,n.fold)
  res[1:3] <- metrics(U,u1,Pi)
  
  ## For proposed
  u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  res[4:6] <- metrics(U,u2,Pi)
  
  return(res)
}

len.stat <- 6
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./primary_results/part2/setting3/","stats_random2.csv"),row.names = FALSE)

## random graph 3
simu <- function(i){
  p <- 10;n <- 200;
  tau.list <- seq(0.1,0.2,0.01)
  gamma.list <- seq(0.05,0.5,0.05)
  n.fold <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(6)
  
  ## Data generation
  generation <- random.generation(p,n,U.test,type)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  W <- generation$W
  V <- W%*%solve(diag(p)- U)
  
  out <- topological_order(V)
  Pi <- out$an_mat
  Phi <- out$in_mat
  Piv <- out$iv_mat
  
  
  ## For direct computing
  u1 <- coef.direct(X,Y,Pi,Phi,tau.list,gamma.list,n.fold)
  res[1:3] <- metrics(U,u1,Pi)
  
  ## For proposed
  u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  res[4:6] <- metrics(U,u2,Pi)
  
  return(res)
}

len.stat <- 6
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./primary_results/part2/setting3/","stats_random3.csv"),row.names = FALSE)

## random graph 4
simu <- function(i){
  p <- 10;n <- 1000;
  tau.list <- seq(0.1,0.2,0.01)
  gamma.list <- seq(0.05,0.5,0.05)
  n.fold <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(6)
  
  ## Data generation
  generation <- random.generation2(p,n,U.test,type)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  W <- generation$W
  V <- W%*%solve(diag(p)- U)
  
  out <- topological_order(V)
  Pi <- out$an_mat
  Phi <- out$in_mat
  Piv <- out$iv_mat
  
  
  ## For direct computing
  u1 <- coef.direct(X,Y,Pi,Phi,tau.list,gamma.list,n.fold)
  res[1:3] <- metrics(U,u1,Pi)
  
  ## For proposed
  u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  res[4:6] <- metrics(U,u2,Pi)
  
  return(res)
}

len.stat <- 6
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./primary_results/part2/setting3/","stats_random4.csv"),row.names = FALSE)

## random graph 5
simu <- function(i){
  p <- 10;n <- 500;
  tau.list <- seq(0.1,0.2,0.01)
  gamma.list <- seq(0.05,0.5,0.05)
  n.fold <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(6)
  
  ## Data generation
  generation <- random.generation2(p,n,U.test,type)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  W <- generation$W
  V <- W%*%solve(diag(p)- U)
  
  out <- topological_order(V)
  Pi <- out$an_mat
  Phi <- out$in_mat
  Piv <- out$iv_mat
  
  
  ## For direct computing
  u1 <- coef.direct(X,Y,Pi,Phi,tau.list,gamma.list,n.fold)
  res[1:3] <- metrics(U,u1,Pi)
  
  ## For proposed
  u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  res[4:6] <- metrics(U,u2,Pi)
  
  return(res)
}

len.stat <- 6
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./primary_results/part2/setting3/","stats_random5.csv"),row.names = FALSE)

## random graph 6
simu <- function(i){
  p <- 10;n <- 200;
  tau.list <- seq(0.1,0.2,0.01)
  gamma.list <- seq(0.05,0.5,0.05)
  n.fold <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(6)
  
  ## Data generation
  generation <- random.generation2(p,n,U.test,type)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  W <- generation$W
  V <- W%*%solve(diag(p)- U)
  
  out <- topological_order(V)
  Pi <- out$an_mat
  Phi <- out$in_mat
  Piv <- out$iv_mat
  
  
  ## For direct computing
  u1 <- coef.direct(X,Y,Pi,Phi,tau.list,gamma.list,n.fold)
  res[1:3] <- metrics(U,u1,Pi)
  
  ## For proposed
  u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  res[4:6] <- metrics(U,u2,Pi)
  
  return(res)
}

len.stat <- 6
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./primary_results/part2/setting3/","stats_random6.csv"),row.names = FALSE)

## hub graph 1
simu <- function(i){
  p <- 11;n <- 1000;
  tau.list <- seq(0.2,0.3,0.01)
  gamma.list <- seq(0.1,1,0.1)
  n.fold <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(6)
  
  ## Data generation
  generation <- hub.generation(p,n)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  W <- generation$W
  V <- W%*%solve(diag(p)- U)
  
  out <- topological_order(V)
  Pi <- out$an_mat
  Phi <- out$in_mat
  Piv <- out$iv_mat
  
  
  ## For direct computing
  u1 <- coef.direct(X,Y,Pi,Phi,tau.list,gamma.list,n.fold)
  res[1:3] <- metrics(U,u1,Pi)
  
  ## For proposed
  u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  res[4:6] <- metrics(U,u2,Pi)
  
  return(res)
}

len.stat <- 6
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);stat <- simu(i);return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./primary_results/part2/setting3/","stats_hub1.csv"),row.names = FALSE)

## hub graph 2
simu <- function(i){
  p <- 11;n <- 500;
  tau.list <- seq(0.2,0.3,0.01)
  gamma.list <- seq(0.1,1,0.1)
  n.fold <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(6)
  
  ## Data generation
  generation <- hub.generation(p,n)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  W <- generation$W
  V <- W%*%solve(diag(p)- U)
  
  out <- topological_order(V)
  Pi <- out$an_mat
  Phi <- out$in_mat
  Piv <- out$iv_mat
  
  
  ## For direct computing
  u1 <- coef.direct(X,Y,Pi,Phi,tau.list,gamma.list,n.fold)
  res[1:3] <- metrics(U,u1,Pi)
  
  ## For proposed
  u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  res[4:6] <- metrics(U,u2,Pi)
  
  return(res)
}

len.stat <- 6
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./primary_results/part2/setting3/","stats_hub2.csv"),row.names = FALSE)

## hub graph 3
simu <- function(i){
  p <- 11;n <- 200;
  tau.list <- seq(0.2,0.3,0.01)
  gamma.list <- seq(0.1,1,0.1)
  n.fold <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(6)
  
  ## Data generation
  generation <- hub.generation(p,n)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  W <- generation$W
  V <- W%*%solve(diag(p)- U)
  
  out <- topological_order(V)
  Pi <- out$an_mat
  Phi <- out$in_mat
  Piv <- out$iv_mat
  
  
  ## For direct computing
  u1 <- coef.direct(X,Y,Pi,Phi,tau.list,gamma.list,n.fold)
  res[1:3] <- metrics(U,u1,Pi)
  
  ## For proposed
  u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  res[4:6] <- metrics(U,u2,Pi)
  
  return(res)
}

len.stat <- 6
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./primary_results/part2/setting3/","stats_hub3.csv"),row.names = FALSE)

## hub graph 4
simu <- function(i){
  p <- 11;n <- 1000;
  tau.list <- seq(0.2,0.3,0.01)
  gamma.list <- seq(0.1,1,0.1)
  n.fold <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(6)
  
  ## Data generation
  generation <- hub.generation2(p,n)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  W <- generation$W
  V <- W%*%solve(diag(p)- U)
  
  out <- topological_order(V)
  Pi <- out$an_mat
  Phi <- out$in_mat
  Piv <- out$iv_mat
  
  
  ## For direct computing
  u1 <- coef.direct(X,Y,Pi,Phi,tau.list,gamma.list,n.fold)
  res[1:3] <- metrics(U,u1,Pi)
  
  ## For proposed
  u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  res[4:6] <- metrics(U,u2,Pi)
  
  return(res)
}

len.stat <- 6
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./primary_results/part2/setting3/","stats_hub4.csv"),row.names = FALSE)

## hub graph 5
simu <- function(i){
  p <- 11;n <- 500;
  tau.list <- seq(0.2,0.3,0.01)
  gamma.list <- seq(0.1,1,0.1)
  n.fold <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(6)
  
  ## Data generation
  generation <- hub.generation2(p,n)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  W <- generation$W
  V <- W%*%solve(diag(p)- U)
  
  out <- topological_order(V)
  Pi <- out$an_mat
  Phi <- out$in_mat
  Piv <- out$iv_mat
  
  
  ## For direct computing
  u1 <- coef.direct(X,Y,Pi,Phi,tau.list,gamma.list,n.fold)
  res[1:3] <- metrics(U,u1,Pi)
  
  ## For proposed
  u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  res[4:6] <- metrics(U,u2,Pi)
  
  return(res)
}

len.stat <- 6
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./primary_results/part2/setting3/","stats_hub5.csv"),row.names = FALSE)

## hub graph 6
simu <- function(i){
  p <- 11;n <- 200;
  tau.list <- seq(0.2,0.3,0.01)
  gamma.list <- seq(0.1,1,0.1)
  n.fold <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(6)
  
  ## Data generation
  generation <- hub.generation2(p,n)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  W <- generation$W
  V <- W%*%solve(diag(p)- U)
  
  out <- topological_order(V)
  Pi <- out$an_mat
  Phi <- out$in_mat
  Piv <- out$iv_mat
  
  
  ## For direct computing
  u1 <- coef.direct(X,Y,Pi,Phi,tau.list,gamma.list,n.fold)
  res[1:3] <- metrics(U,u1,Pi)
  
  ## For proposed
  u2 <- coef.proposed(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  res[4:6] <- metrics(U,u2,Pi)
  
  return(res)
}

len.stat <- 6
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./primary_results/part2/setting3/","stats_hub6.csv"),row.names = FALSE)

stopCluster(cl)