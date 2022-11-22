library(doParallel)
n_cores <- detectCores() - 1
registerDoParallel(cores=n_cores)  
cl <- makeCluster(n_cores, type="FORK")

clusterEvalQ(cl,{source("../mbtlp.R",chdir = TRUE)
  source("../intdagdiscovery.R",chdir = TRUE)
  source("../intdagcoef.R",chdir = TRUE)
  source("../intdaginfer.r",chdir = TRUE)
  library(mvtnorm)
  library(MASS)
  library(clusterGeneration)
  library(mnormt)})

simu.single.2lr <- function(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list,max.it=100000,tol=1e-7){
  
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

random.generation <- function(p,n,U.test,type){
  if (p%%2!=0)
    stop("wrong p is given")
  q <- 2.5*p
  U <- matrix(sample(x=c(1,0),size=p*p,replace=TRUE,prob=c(1/p,1-1/p)),p,p)
  for (i in 1:p){
    for (j in 1:i){
      U[i,j] <- 0
    }
  }
  W1 <- W2 <- diag(1,p,p); W3 <- matrix(0,0.5*p,p); 
  for (i in 1:(p/2)){
    W3[i,2*i-1] <- 1
    W3[i,2*i] <- 1
  }
  W <- rbind(W1,W2,W3)
  if (type==0){
    U <- (1-U.test)*U
  }else{
    U <- 1*((U+U.test)>0)
  }
  U <- U*matrix(runif(p*p,min=0.8,max=1.2)*sample(c(-1,1),p*p,replace = TRUE),p,p)
  W <- W*matrix(runif(q*p,min=0.8,max=1.2)*sample(c(-1,1),q*p,replace = TRUE),q,p)
  vari <- diag(runif(p,min=0.4,max=0.6))
  sigma <- vari%*%diag(p)%*%vari
  e <- mvrnorm(n,rep(0,p),sigma)
  r <- p/2
  V <- matrix(0,r,p)
  for (i in 1:r){
    V[i,2*i-1] <- 1
    V[i,2*i] <- 1
  }
  V <- V*matrix(runif(r*p,min=0.4,max=0.6)*sample(c(-1,1),r*p,replace = TRUE),r,p)
  Z <- matrix(rnorm(r*n),ncol=r)
  X <- matrix(rnorm(q*n),ncol=q); Y <- matrix(0,n,p); Y <- (X%*%W+e+Z%*%V)%*%solve(diag(1,p)-U);
  Sigma <- sigma+t(V)%*%V
  return(list(X=X,Y=Y,U=U,W=W,Sigma=Sigma))
}

random.generation2 <- function(p,n,U.test,type){
  if (p%%2!=0)
    stop("wrong p is given")
  q <- 2.5*p
  U <- matrix(sample(x=c(1,0),size=p*p,replace=TRUE,prob=c(1/p,1-1/p)),p,p)
  for (i in 1:p){
    for (j in 1:i){
      U[i,j] <- 0
    }
  }
  W1 <- W2 <- diag(1,p,p); W3 <- matrix(0,0.5*p,p); 
  for (i in 1:(p/2)){
    W3[i,2*i-1] <- 1
    W3[i,2*i] <- 1
  }
  W <- rbind(W1,W2,W3)
  if (type==0){
    U <- (1-U.test)*U
  }else{
    U <- 1*((U+U.test)>0)
  }
  U <- U*matrix(runif(p*p,min=0.8,max=1.2)*sample(c(-1,1),p*p,replace = TRUE),p,p)
  W <- W*matrix(runif(q*p,min=0.8,max=1.2)*sample(c(-1,1),q*p,replace = TRUE),q,p)
  vari <- diag(runif(p,min=0.4,max=0.6))
  sigma <- vari%*%diag(p)%*%vari
  e <- mvrnorm(n,rep(0,p),sigma)
  r <- p/2
  V <- matrix(0,r,p)
  for (i in 1:r){
    V[i,2*i-1] <- 1
    V[i,2*i] <- 1
  }
  V <- V*matrix(runif(r*p,min=0.4,max=0.6)*sample(c(-1,1),r*p,replace = TRUE),r,p)
  Z <- matrix(rnorm(r*n),ncol=r)
  X <- matrix(sample(c(-1,1),q*n,replace = TRUE),ncol=q); Y <- matrix(0,n,p); Y <- (X%*%W+e+Z%*%V)%*%solve(diag(1,p)-U);
  Sigma <- sigma+t(V)%*%V
  return(list(X=X,Y=Y,U=U,W=W,Sigma=Sigma))
}


simu.random <- function(times,p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list,type){
  stats <- matrix(0,times,length(U.test.list))
  U.test.bind <- matrix(0,p,p)
  for (k in 1:length(U.test.list)){
    U.test.bind <- U.test.bind + U.test.list[[k]]
  }
  U.test.bind <- 1*(U.test.bind >0)
  for (i in 1:times){
    generation <- random.generation(p,n,U.test.bind,type)
    X <- generation$X
    Y <- generation$Y
    stats[i,] <- simu.single.2lr(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list)
  }
  return(stats)
}

simu.random2 <- function(times,p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list,type){
  stats <- matrix(0,times,length(U.test.list))
  U.test.bind <- matrix(0,p,p)
  for (k in 1:length(U.test.list)){
    U.test.bind <- U.test.bind + U.test.list[[k]]
  }
  U.test.bind <- 1*(U.test.bind >0)
  for (i in 1:times){
    generation <- random.generation2(p,n,U.test.bind,type)
    X <- generation$X
    Y <- generation$Y
    stats[i,] <- simu.single.2lr(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list)
  }
  return(stats)
}

clusterExport(cl,c("simu.single.2lr","is.acyclic","random.generation","random.generation2","simu.random","simu.random2"))

## random graph 1
simu <- function(i){
  times <- 1;p <- 10;n <- 1200;
  tau.list1 <- seq(0.1,0.2,0.01)
  gamma.list1 <- seq(0.05,0.5,0.05)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  tau.list3 <- seq(0.05,0.1,0.01)
  gamma.list3 <- seq(0.005,0.1,0.005)
  n.fold1 <- n.fold2 <- n.fold3 <- 5
  U.test.list1 <- vector("list",3)
  U.test.list2 <- vector("list",3)
  U.test1 <- matrix(0,p,p)
  U.test1[1,3] <- 1;
  U.test2 <- matrix(0,p,p)
  U.test2[1,3] <- 1; U.test2[2,4] <- 1; U.test2[5,7] <- 1;
  U.test3 <- matrix(0,p,p)
  U.test3[1,3] <- 1; U.test3[2,4] <- 1; U.test3[5,7] <- 1; U.test3[6,8] <- 1; U.test3[9,10] <- 1;
  U.test4 <- matrix(0,p,p)
  U.test4[1,3] <- 1;
  U.test5 <- matrix(0,p,p)
  U.test5[1,3] <- 1; U.test5[2,4] <- 1; U.test5[5,7] <- 1;
  U.test6 <- matrix(0,p,p)
  U.test6[1,3] <- 1; U.test6[2,4] <- 1; U.test6[5,7] <- 1; U.test6[6,8] <- 1; U.test6[9,10] <- 1; 
  U.test.list1[[1]] <- U.test1; U.test.list1[[2]] <- U.test2; U.test.list1[[3]] <- U.test3;
  U.test.list2[[1]] <- U.test4; U.test.list2[[2]] <- U.test5; U.test.list2[[3]] <- U.test6;
  stats.random <- numeric(6)
  stats.random[1:3] <- simu.random(times,p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list1,type=0)
  stats.random[4:6] <- simu.random(times,p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list2,type=1)
  return(stats.random)
}

len.test <- 6
clusterExport(cl,c("simu","len.test"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.test));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./part3/setting2/","stats_random1.csv"),row.names = FALSE)

## random graph 2
simu <- function(i){
  times <- 1;p <- 10;n <- 700;
  tau.list1 <- seq(0.1,0.2,0.01)
  gamma.list1 <- seq(0.05,0.5,0.05)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  tau.list3 <- seq(0.05,0.1,0.01)
  gamma.list3 <- seq(0.005,0.1,0.005)
  n.fold1 <- n.fold2 <- n.fold3 <- 5
  U.test.list1 <- vector("list",3)
  U.test.list2 <- vector("list",3)
  U.test1 <- matrix(0,p,p)
  U.test1[1,3] <- 1;
  U.test2 <- matrix(0,p,p)
  U.test2[1,3] <- 1; U.test2[2,4] <- 1; U.test2[5,7] <- 1;
  U.test3 <- matrix(0,p,p)
  U.test3[1,3] <- 1; U.test3[2,4] <- 1; U.test3[5,7] <- 1; U.test3[6,8] <- 1; U.test3[9,10] <- 1;
  U.test4 <- matrix(0,p,p)
  U.test4[1,3] <- 1;
  U.test5 <- matrix(0,p,p)
  U.test5[1,3] <- 1; U.test5[2,4] <- 1; U.test5[5,7] <- 1;
  U.test6 <- matrix(0,p,p)
  U.test6[1,3] <- 1; U.test6[2,4] <- 1; U.test6[5,7] <- 1; U.test6[6,8] <- 1; U.test6[9,10] <- 1; 
  U.test.list1[[1]] <- U.test1; U.test.list1[[2]] <- U.test2; U.test.list1[[3]] <- U.test3;
  U.test.list2[[1]] <- U.test4; U.test.list2[[2]] <- U.test5; U.test.list2[[3]] <- U.test6;
  stats.random <- numeric(6)
  stats.random[1:3] <- simu.random(times,p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list1,type=0)
  stats.random[4:6] <- simu.random(times,p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list2,type=1)
  return(stats.random)
}

len.test <- 6
clusterExport(cl,c("simu","len.test"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.test));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./part3/setting2/","stats_random2.csv"),row.names = FALSE)

## random graph 3
simu <- function(i){
  times <- 1;p <- 10;n <- 400;
  tau.list1 <- seq(0.1,0.2,0.01)
  gamma.list1 <- seq(0.05,0.5,0.05)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  tau.list3 <- seq(0.05,0.1,0.01)
  gamma.list3 <- seq(0.005,0.1,0.005)
  n.fold1 <- n.fold2 <- n.fold3 <- 5
  U.test.list1 <- vector("list",3)
  U.test.list2 <- vector("list",3)
  U.test1 <- matrix(0,p,p)
  U.test1[1,3] <- 1;
  U.test2 <- matrix(0,p,p)
  U.test2[1,3] <- 1; U.test2[2,4] <- 1; U.test2[5,7] <- 1;
  U.test3 <- matrix(0,p,p)
  U.test3[1,3] <- 1; U.test3[2,4] <- 1; U.test3[5,7] <- 1; U.test3[6,8] <- 1; U.test3[9,10] <- 1;
  U.test4 <- matrix(0,p,p)
  U.test4[1,3] <- 1;
  U.test5 <- matrix(0,p,p)
  U.test5[1,3] <- 1; U.test5[2,4] <- 1; U.test5[5,7] <- 1;
  U.test6 <- matrix(0,p,p)
  U.test6[1,3] <- 1; U.test6[2,4] <- 1; U.test6[5,7] <- 1; U.test6[6,8] <- 1; U.test6[9,10] <- 1; 
  U.test.list1[[1]] <- U.test1; U.test.list1[[2]] <- U.test2; U.test.list1[[3]] <- U.test3;
  U.test.list2[[1]] <- U.test4; U.test.list2[[2]] <- U.test5; U.test.list2[[3]] <- U.test6;
  stats.random <- numeric(6)
  stats.random[1:3] <- simu.random(times,p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list1,type=0)
  stats.random[4:6] <- simu.random(times,p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list2,type=1)
  return(stats.random)
}

len.test <- 6
clusterExport(cl,c("simu","len.test"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.test));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./part3/setting2/","stats_random3.csv"),row.names = FALSE)

## random graph 4
simu <- function(i){
  times <- 1;p <- 10;n <- 1200;
  tau.list1 <- seq(0.1,0.2,0.01)
  gamma.list1 <- seq(0.05,0.5,0.05)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  tau.list3 <- seq(0.05,0.1,0.01)
  gamma.list3 <- seq(0.005,0.1,0.005)
  n.fold1 <- n.fold2 <- n.fold3 <- 5
  U.test.list1 <- vector("list",3)
  U.test.list2 <- vector("list",3)
  U.test1 <- matrix(0,p,p)
  U.test1[1,3] <- 1;
  U.test2 <- matrix(0,p,p)
  U.test2[1,3] <- 1; U.test2[2,4] <- 1; U.test2[5,7] <- 1;
  U.test3 <- matrix(0,p,p)
  U.test3[1,3] <- 1; U.test3[2,4] <- 1; U.test3[5,7] <- 1; U.test3[6,8] <- 1; U.test3[9,10] <- 1;
  U.test4 <- matrix(0,p,p)
  U.test4[1,3] <- 1;
  U.test5 <- matrix(0,p,p)
  U.test5[1,3] <- 1; U.test5[2,4] <- 1; U.test5[5,7] <- 1;
  U.test6 <- matrix(0,p,p)
  U.test6[1,3] <- 1; U.test6[2,4] <- 1; U.test6[5,7] <- 1; U.test6[6,8] <- 1; U.test6[9,10] <- 1; 
  U.test.list1[[1]] <- U.test1; U.test.list1[[2]] <- U.test2; U.test.list1[[3]] <- U.test3;
  U.test.list2[[1]] <- U.test4; U.test.list2[[2]] <- U.test5; U.test.list2[[3]] <- U.test6;
  stats.random <- numeric(6)
  stats.random[1:3] <- simu.random2(times,p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list1,type=0)
  stats.random[4:6] <- simu.random2(times,p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list2,type=1)
  return(stats.random)
}

len.test <- 6
clusterExport(cl,c("simu","len.test"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.test));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./part3/setting2/","stats_random4.csv"),row.names = FALSE)


## random graph 5
simu <- function(i){
  times <- 1;p <- 10;n <- 700;
  tau.list1 <- seq(0.1,0.2,0.01)
  gamma.list1 <- seq(0.05,0.5,0.05)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  tau.list3 <- seq(0.05,0.1,0.01)
  gamma.list3 <- seq(0.005,0.1,0.005)
  n.fold1 <- n.fold2 <- n.fold3 <- 5
  U.test.list1 <- vector("list",3)
  U.test.list2 <- vector("list",3)
  U.test1 <- matrix(0,p,p)
  U.test1[1,3] <- 1;
  U.test2 <- matrix(0,p,p)
  U.test2[1,3] <- 1; U.test2[2,4] <- 1; U.test2[5,7] <- 1;
  U.test3 <- matrix(0,p,p)
  U.test3[1,3] <- 1; U.test3[2,4] <- 1; U.test3[5,7] <- 1; U.test3[6,8] <- 1; U.test3[9,10] <- 1;
  U.test4 <- matrix(0,p,p)
  U.test4[1,3] <- 1;
  U.test5 <- matrix(0,p,p)
  U.test5[1,3] <- 1; U.test5[2,4] <- 1; U.test5[5,7] <- 1;
  U.test6 <- matrix(0,p,p)
  U.test6[1,3] <- 1; U.test6[2,4] <- 1; U.test6[5,7] <- 1; U.test6[6,8] <- 1; U.test6[9,10] <- 1; 
  U.test.list1[[1]] <- U.test1; U.test.list1[[2]] <- U.test2; U.test.list1[[3]] <- U.test3;
  U.test.list2[[1]] <- U.test4; U.test.list2[[2]] <- U.test5; U.test.list2[[3]] <- U.test6;
  stats.random <- numeric(6)
  stats.random[1:3] <- simu.random2(times,p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list1,type=0)
  stats.random[4:6] <- simu.random2(times,p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list2,type=1)
  return(stats.random)
}

len.test <- 6
clusterExport(cl,c("simu","len.test"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.test));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./part3/setting2/","stats_random5.csv"),row.names = FALSE)

## random graph 6
simu <- function(i){
  times <- 1;p <- 10;n <- 400;
  tau.list1 <- seq(0.1,0.2,0.01)
  gamma.list1 <- seq(0.05,0.5,0.05)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  tau.list3 <- seq(0.05,0.1,0.01)
  gamma.list3 <- seq(0.005,0.1,0.005)
  n.fold1 <- n.fold2 <- n.fold3 <- 5
  U.test.list1 <- vector("list",3)
  U.test.list2 <- vector("list",3)
  U.test1 <- matrix(0,p,p)
  U.test1[1,3] <- 1;
  U.test2 <- matrix(0,p,p)
  U.test2[1,3] <- 1; U.test2[2,4] <- 1; U.test2[5,7] <- 1;
  U.test3 <- matrix(0,p,p)
  U.test3[1,3] <- 1; U.test3[2,4] <- 1; U.test3[5,7] <- 1; U.test3[6,8] <- 1; U.test3[9,10] <- 1;
  U.test4 <- matrix(0,p,p)
  U.test4[1,3] <- 1;
  U.test5 <- matrix(0,p,p)
  U.test5[1,3] <- 1; U.test5[2,4] <- 1; U.test5[5,7] <- 1;
  U.test6 <- matrix(0,p,p)
  U.test6[1,3] <- 1; U.test6[2,4] <- 1; U.test6[5,7] <- 1; U.test6[6,8] <- 1; U.test6[9,10] <- 1; 
  U.test.list1[[1]] <- U.test1; U.test.list1[[2]] <- U.test2; U.test.list1[[3]] <- U.test3;
  U.test.list2[[1]] <- U.test4; U.test.list2[[2]] <- U.test5; U.test.list2[[3]] <- U.test6;
  stats.random <- numeric(6)
  stats.random[1:3] <- simu.random2(times,p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list1,type=0)
  stats.random[4:6] <- simu.random2(times,p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list2,type=1)
  return(stats.random)
}

len.test <- 6
clusterExport(cl,c("simu","len.test"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.test));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=6)
write.csv(stats_random,file.path("./part3/setting2/","stats_random6.csv"),row.names = FALSE)

stopCluster(cl)