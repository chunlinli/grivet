library(doParallel)
n_cores <- detectCores() - 1
registerDoParallel(cores=n_cores)  
cl <- makeCluster(n_cores, type="FORK")

clusterEvalQ(cl,{source("comparison1_setting2.R",chdir = TRUE)
  library(mvtnorm)
  library(MASS)})

## random graph 1
simu <- function(i){
  p <- 10;n <- 1000;
  tau.list1 <- seq(0.1,0.2,0.01)
  gamma.list1 <- seq(0.05,0.5,0.05)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  n.fold1 <- n.fold2 <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(12)
  
  ## Data generation
  generation <- random.generation(p,n,U.test,type)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  U <- 1*(U!=0)
  y <- cbind(Y,X)
  
  
  ## For RFCI
  out1 <- causal_discovery_RFCI(y)
  u1 <- out1$u
  u1 <- u1[1:p,1:p]
  
  u1 <- PAG2best(U,u1)
  res[1:4] <- metrics(U,u1)
  
  ## For Lrps-ges
  out2 <- causal_discovery_LRpS_GES(y)
  u2 <- out2$u
  u2 <- u2[1:p,1:p]
  res[5:8] <- metrics(U,u2)
  
  ## For proposed
  u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
  res[9:12] <- metrics(U,u3)
  
  return(res)
}

len.stat <- 12
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);stat<-simu(i);return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=12)
write.csv(stats_random,file.path("./part1/setting2/","stats_random1.csv"),row.names = FALSE)

## random graph 2
simu <- function(i){
  p <- 10;n <- 500;
  tau.list1 <- seq(0.1,0.2,0.01)
  gamma.list1 <- seq(0.05,0.5,0.05)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  n.fold1 <- n.fold2 <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(12)
  
  ## Data generation
  generation <- random.generation(p,n,U.test,type)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  U <- 1*(U!=0)
  y <- cbind(Y,X)
  
  
  ## For RFCI
  out1 <- causal_discovery_RFCI(y)
  u1 <- out1$u
  u1 <- u1[1:p,1:p]
  
  u1 <- PAG2best(U,u1)
  res[1:4] <- metrics(U,u1)
  
  ## For Lrps-ges
  out2 <- causal_discovery_LRpS_GES(y)
  u2 <- out2$u
  u2 <- u2[1:p,1:p]
  res[5:8] <- metrics(U,u2)
  
  ## For proposed
  u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
  res[9:12] <- metrics(U,u3)
  
  return(res)
}

len.stat <- 12
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=12)
write.csv(stats_random,file.path("./part1/setting2/","stats_random2.csv"),row.names = FALSE)


## random graph 3
simu <- function(i){
  p <- 10;n <- 200;
  tau.list1 <- seq(0.1,0.2,0.01)
  gamma.list1 <- seq(0.05,0.5,0.05)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  n.fold1 <- n.fold2 <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(12)
  
  ## Data generation
  generation <- random.generation(p,n,U.test,type)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  U <- 1*(U!=0)
  y <- cbind(Y,X)
  
  
  ## For RFCI
  out1 <- causal_discovery_RFCI(y)
  u1 <- out1$u
  u1 <- u1[1:p,1:p]
  
  u1 <- PAG2best(U,u1)
  res[1:4] <- metrics(U,u1)
  
  ## For Lrps-ges
  out2 <- causal_discovery_LRpS_GES(y)
  u2 <- out2$u
  u2 <- u2[1:p,1:p]
  res[5:8] <- metrics(U,u2)
  
  ## For proposed
  u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
  res[9:12] <- metrics(U,u3)
  
  return(res)
}

len.stat <- 12
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=12)
write.csv(stats_random,file.path("./part1/setting2/","stats_random3.csv"),row.names = FALSE)

## random graph 4
simu <- function(i){
  p <- 10;n <- 1000;
  tau.list1 <- seq(0.1,0.2,0.01)
  gamma.list1 <- seq(0.05,0.5,0.05)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  n.fold1 <- n.fold2 <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(12)
  
  ## Data generation
  generation <- random.generation2(p,n,U.test,type)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  U <- 1*(U!=0)
  y <- cbind(Y,X)
  
  
  ## For RFCI
  out1 <- causal_discovery_RFCI(y)
  u1 <- out1$u
  u1 <- u1[1:p,1:p]
  
  u1 <- PAG2best(U,u1)
  res[1:4] <- metrics(U,u1)
  
  ## For Lrps-ges
  out2 <- causal_discovery_LRpS_GES(y)
  u2 <- out2$u
  u2 <- u2[1:p,1:p]
  res[5:8] <- metrics(U,u2)
  
  ## For proposed
  u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
  res[9:12] <- metrics(U,u3)
  
  return(res)
}

len.stat <- 12
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=12)
write.csv(stats_random,file.path("./part1/setting2/","stats_random4.csv"),row.names = FALSE)

## random graph 5
simu <- function(i){
  p <- 10;n <- 500;
  tau.list1 <- seq(0.1,0.2,0.01)
  gamma.list1 <- seq(0.05,0.5,0.05)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  n.fold1 <- n.fold2 <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(12)
  
  ## Data generation
  generation <- random.generation2(p,n,U.test,type)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  U <- 1*(U!=0)
  y <- cbind(Y,X)
  
  
  ## For RFCI
  out1 <- causal_discovery_RFCI(y)
  u1 <- out1$u
  u1 <- u1[1:p,1:p]
  
  u1 <- PAG2best(U,u1)
  res[1:4] <- metrics(U,u1)
  
  ## For Lrps-ges
  out2 <- causal_discovery_LRpS_GES(y)
  u2 <- out2$u
  u2 <- u2[1:p,1:p]
  res[5:8] <- metrics(U,u2)
  
  ## For proposed
  u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
  res[9:12] <- metrics(U,u3)
  
  return(res)
}

len.stat <- 12
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=12)
write.csv(stats_random,file.path("./part1/setting2/","stats_random5.csv"),row.names = FALSE)


## random graph 6
simu <- function(i){
  p <- 10;n <- 200;
  tau.list1 <- seq(0.1,0.2,0.01)
  gamma.list1 <- seq(0.05,0.5,0.05)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  n.fold1 <- n.fold2 <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(12)
  
  ## Data generation
  generation <- random.generation2(p,n,U.test,type)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  U <- 1*(U!=0)
  y <- cbind(Y,X)
  
  
  ## For RFCI
  out1 <- causal_discovery_RFCI(y)
  u1 <- out1$u
  u1 <- u1[1:p,1:p]
  
  u1 <- PAG2best(U,u1)
  res[1:4] <- metrics(U,u1)
  
  ## For Lrps-ges
  out2 <- causal_discovery_LRpS_GES(y)
  u2 <- out2$u
  u2 <- u2[1:p,1:p]
  res[5:8] <- metrics(U,u2)
  
  ## For proposed
  u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
  res[9:12] <- metrics(U,u3)
  
  return(res)
}

len.stat <- 12
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=12)
write.csv(stats_random,file.path("./part1/setting2/","stats_random6.csv"),row.names = FALSE)

## hub graph 1
simu <- function(i){
  p <- 11;n <- 1000;
  tau.list1 <- seq(0.2,0.3,0.01)
  gamma.list1 <- seq(0.1,1,0.1)
  tau.list2 <- seq(0.2,0.3,0.01)
  gamma.list2 <- seq(0.1,1,0.1)
  n.fold1 <- n.fold2 <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(12)
  
  ## Data generation
  generation <- hub.generation(p,n)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  U <- 1*(U!=0)
  y <- cbind(Y,X)
  
  
  ## For RFCI
  out1 <- causal_discovery_RFCI(y)
  u1 <- out1$u
  u1 <- u1[1:p,1:p]
  
  u1 <- PAG2best(U,u1)
  res[1:4] <- metrics(U,u1)
  
  ## For Lrps-ges
  out2 <- causal_discovery_LRpS_GES(y)
  u2 <- out2$u
  u2 <- u2[1:p,1:p]
  res[5:8] <- metrics(U,u2)
  
  ## For proposed
  u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
  res[9:12] <- metrics(U,u3)
  
  return(res)
}

len.stat <- 12
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=12)
write.csv(stats_random,file.path("./part1/setting2/","stats_hub1.csv"),row.names = FALSE)

## hub graph 2
simu <- function(i){
  p <- 11;n <- 500;
  tau.list1 <- seq(0.2,0.3,0.01)
  gamma.list1 <- seq(0.1,1,0.1)
  tau.list2 <- seq(0.2,0.3,0.01)
  gamma.list2 <- seq(0.1,1,0.1)
  n.fold1 <- n.fold2 <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(12)
  
  ## Data generation
  generation <- hub.generation(p,n)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  U <- 1*(U!=0)
  y <- cbind(Y,X)
  
  
  ## For RFCI
  out1 <- causal_discovery_RFCI(y)
  u1 <- out1$u
  u1 <- u1[1:p,1:p]
  
  u1 <- PAG2best(U,u1)
  res[1:4] <- metrics(U,u1)
  
  ## For Lrps-ges
  out2 <- causal_discovery_LRpS_GES(y)
  u2 <- out2$u
  u2 <- u2[1:p,1:p]
  res[5:8] <- metrics(U,u2)
  
  ## For proposed
  u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
  res[9:12] <- metrics(U,u3)
  
  return(res)
}

len.stat <- 12
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=12)
write.csv(stats_random,file.path("./part1/setting2/","stats_hub2.csv"),row.names = FALSE)


## hub graph 3
simu <- function(i){
  p <- 11;n <- 200;
  tau.list1 <- seq(0.2,0.3,0.01)
  gamma.list1 <- seq(0.1,1,0.1)
  tau.list2 <- seq(0.2,0.3,0.01)
  gamma.list2 <- seq(0.1,1,0.1)
  n.fold1 <- n.fold2 <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(12)
  
  ## Data generation
  generation <- hub.generation(p,n)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  U <- 1*(U!=0)
  y <- cbind(Y,X)
  
  
  ## For RFCI
  out1 <- causal_discovery_RFCI(y)
  u1 <- out1$u
  u1 <- u1[1:p,1:p]
  
  u1 <- PAG2best(U,u1)
  res[1:4] <- metrics(U,u1)
  
  ## For Lrps-ges
  out2 <- causal_discovery_LRpS_GES(y)
  u2 <- out2$u
  u2 <- u2[1:p,1:p]
  res[5:8] <- metrics(U,u2)
  
  ## For proposed
  u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
  res[9:12] <- metrics(U,u3)
  
  return(res)
}

len.stat <- 12
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=12)
write.csv(stats_random,file.path("./part1/setting2/","stats_hub3.csv"),row.names = FALSE)

## hub graph 4
simu <- function(i){
  p <- 11;n <- 1000;
  tau.list1 <- seq(0.2,0.3,0.01)
  gamma.list1 <- seq(0.1,1,0.1)
  tau.list2 <- seq(0.2,0.3,0.01)
  gamma.list2 <- seq(0.1,1,0.1)
  n.fold1 <- n.fold2 <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(12)
  
  ## Data generation
  generation <- hub.generation2(p,n)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  U <- 1*(U!=0)
  y <- cbind(Y,X)
  
  
  ## For RFCI
  out1 <- causal_discovery_RFCI(y)
  u1 <- out1$u
  u1 <- u1[1:p,1:p]
  
  u1 <- PAG2best(U,u1)
  res[1:4] <- metrics(U,u1)
  
  ## For Lrps-ges
  out2 <- causal_discovery_LRpS_GES(y)
  u2 <- out2$u
  u2 <- u2[1:p,1:p]
  res[5:8] <- metrics(U,u2)
  
  ## For proposed
  u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
  res[9:12] <- metrics(U,u3)
  
  return(res)
}

len.stat <- 12
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=12)
write.csv(stats_random,file.path("./part1/setting2/","stats_hub4.csv"),row.names = FALSE)

## hub graph 5
simu <- function(i){
  p <- 11;n <- 500;
  tau.list1 <- seq(0.2,0.3,0.01)
  gamma.list1 <- seq(0.1,1,0.1)
  tau.list2 <- seq(0.2,0.3,0.01)
  gamma.list2 <- seq(0.1,1,0.1)
  n.fold1 <- n.fold2 <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(12)
  
  ## Data generation
  generation <- hub.generation2(p,n)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  U <- 1*(U!=0)
  y <- cbind(Y,X)
  
  
  ## For RFCI
  out1 <- causal_discovery_RFCI(y)
  u1 <- out1$u
  u1 <- u1[1:p,1:p]
  
  u1 <- PAG2best(U,u1)
  res[1:4] <- metrics(U,u1)
  
  ## For Lrps-ges
  out2 <- causal_discovery_LRpS_GES(y)
  u2 <- out2$u
  u2 <- u2[1:p,1:p]
  res[5:8] <- metrics(U,u2)
  
  ## For proposed
  u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
  res[9:12] <- metrics(U,u3)
  
  return(res)
}

len.stat <- 12
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=12)
write.csv(stats_random,file.path("./part1/setting2/","stats_hub5.csv"),row.names = FALSE)


## hub graph 6
simu <- function(i){
  p <- 11;n <- 200;
  tau.list1 <- seq(0.2,0.3,0.01)
  gamma.list1 <- seq(0.1,1,0.1)
  tau.list2 <- seq(0.2,0.3,0.01)
  gamma.list2 <- seq(0.1,1,0.1)
  n.fold1 <- n.fold2 <- 5
  U.test <- matrix(0,p,p)
  type <- 0
  res <- numeric(12)
  
  ## Data generation
  generation <- hub.generation2(p,n)
  X <- generation$X
  Y <- generation$Y
  U <- generation$U
  U <- 1*(U!=0)
  y <- cbind(Y,X)
  
  
  ## For RFCI
  out1 <- causal_discovery_RFCI(y)
  u1 <- out1$u
  u1 <- u1[1:p,1:p]
  
  u1 <- PAG2best(U,u1)
  res[1:4] <- metrics(U,u1)
  
  ## For Lrps-ges
  out2 <- causal_discovery_LRpS_GES(y)
  u2 <- out2$u
  u2 <- u2[1:p,1:p]
  res[5:8] <- metrics(U,u2)
  
  ## For proposed
  u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
  res[9:12] <- metrics(U,u3)
  
  return(res)
}

len.stat <- 12
clusterExport(cl,c("simu","len.stat"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=12)
write.csv(stats_random,file.path("./part1/setting2/","stats_hub6.csv"),row.names = FALSE)

stopCluster(cl)