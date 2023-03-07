library(doParallel)
n_cores <- detectCores() - 1
registerDoParallel(cores=n_cores)  
cl <- makeCluster(n_cores, type="FORK")

clusterEvalQ(cl,{library(grivet)
  library(mvtnorm)
  library(MASS)
  library(clusterGeneration)
  library(mnormt)})

test.generation <- function(p,len){
  idx <- sample(1:(0.5*p*(p-1)),len)
  iter <- 0; k <- 1; U.test.list <- vector("list",len)
  for(i in 1:(p-1)){
    for(j in (i+1):(p)){
      iter <- iter + 1
      if (iter %in% idx){
        U <- matrix(0,p,p)
        U[i,j] <- 1
        U.test.list[[k]] <- U
        k <- k+1
      }
    }
  }
  return(U.test.list)
}

simu.single.2lr <- function(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list,max.it=100000,tol=1e-7){
  
  n <- nrow(X)
  
  
  # stage 1: discovery
  result.stage1_1 <- cv.intdag.pmle.diff.aic(X,Y,tau.list1,gamma.list1,n.fold1)
  V.es <- result.stage1_1$V
  result.stage1_2 <- topological_order(V.es)
  Pi_es <- result.stage1_2$an_mat
  Phi_es <- result.stage1_2$in_mat
  Piv_es <- result.stage1_2$iv_mat
  
  ## Obtain covariance matrix estimation with respect to Pi_es
  result.stage2 <- cv.intdag.coe(X,Y,Pi_es,Phi_es,Piv_es,tau.list2,gamma.list2,n.fold2)
  Sigma_es <- result.stage2$Sigma
  U_es <- result.stage2$U 
  W_es <- result.stage2$W
  
  
  ## Compute 2lr using the sub-graph
  statistic <- numeric(length(U.test.list))
  for(k in 1:length(U.test.list)){
    U.test <- U.test.list[[k]]
    graph.local <- intdag.localization.ver2(X,Y,U.test,U_es,W_es,Pi_es,Phi_es,Sigma_es)
    U.test.new <- graph.local$U.test.new
    Pi_es.new <- graph.local$Pi.new
    Phi_es.new <- graph.local$Phi.new
    U1.new <- 1*((Pi_es.new + U.test.new)>0)
    if (is.acyclic(U1.new)==0){
      statistic = 0
      return(statistic)
    }else{
      U0.new <- U1.new - U.test.new
    }
    U.hat.new <- graph.local$U.hat.new
    W.hat.new <- graph.local$W.hat.new
    X.new <- graph.local$X.new
    Y.new <- graph.local$Y.new
    Z.new <- Y.new - Y.new%*%U.hat.new - X.new%*%W.hat.new
    MB.support.new <- cv.MB_Union(Z.new,tau.list3,gamma.list3,n.fold3)$S
    Sigma_es.new <- graph.local$Sigma.new
    S_es.new <- precision_refit(Sigma_es.new,MB.support.new,max.it,tol)
    try(result.stage3 <- intdag.2lr(X.new,Y.new,U0.new,Phi_es.new,U1.new,Phi_es.new,S_es.new)$statistic,result.stage3 <- -1)
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

random.generation <- function(p,n,U){
  if (p%%10!=0)
    stop("wrong p is given")
  q <- 2.5*p
  
  W1 <- W2 <- diag(1,p,p); W3 <- matrix(0,0.5*p,p); 
  for (i in 1:(p/2)){
    W3[i,2*i-1] <- 1
    W3[i,2*i] <- 1
  }
  W <- rbind(W1,W2,W3)
  
  U <- U*matrix(rep(1,p*p),p,p)
  W <- W*matrix(rep(1,q*p),q,p)
  vari <- diag(runif(p,min=0.4,max=0.6))
  sigma <- vari%*%diag(p)%*%vari
  e <- mvrnorm(n,rep(0,p),sigma)
  r <- p/10
  V <- matrix(0,r,p)
  for (i in 1:r){
    V[i,(10*(i-1)+1):(10*i)] <- 1
  }
  V <- V*matrix(runif(r*p,min=0.4,max=0.6)*sample(c(-1,1),r*p,replace = TRUE),r,p)
  Z <- matrix(rnorm(r*n),ncol=r)
  X <- matrix(rnorm(q*n),ncol=q); Y <- matrix(0,n,p); Y <- (X%*%W+e+Z%*%V)%*%solve(diag(1,p)-U);
  Sigma <- sigma+t(V)%*%V
  return(list(X=X,Y=Y,U=U,W=W,Sigma=Sigma))
}

random.generation2 <- function(p,n,U){
  if (p%%10!=0)
    stop("wrong p is given")
  q <- 2.5*p
  
  W1 <- W2 <- diag(1,p,p); W3 <- matrix(0,0.5*p,p); 
  for (i in 1:(p/2)){
    W3[i,2*i-1] <- 1
    W3[i,2*i] <- 1
  }
  W <- rbind(W1,W2,W3)
  
  U <- U*matrix(rep(1,p*p),p,p)
  W <- W*matrix(rep(1,q*p),q,p)
  vari <- diag(runif(p,min=0.4,max=0.6))
  sigma <- vari%*%diag(p)%*%vari
  e <- mvrnorm(n,rep(0,p),sigma)
  r <- p/10
  V <- matrix(0,r,p)
  for (i in 1:r){
    V[i,(10*(i-1)+1):(10*i)] <- 1
  }
  V <- V*matrix(runif(r*p,min=0.4,max=0.6)*sample(c(-1,1),r*p,replace = TRUE),r,p)
  Z <- matrix(rnorm(r*n),ncol=r)
  X <- matrix(sample(c(-1,1),q*n,replace = TRUE),ncol=q); Y <- matrix(0,n,p); Y <- (X%*%W+e+Z%*%V)%*%solve(diag(1,p)-U);
  Sigma <- sigma+t(V)%*%V
  return(list(X=X,Y=Y,U=U,W=W,Sigma=Sigma))
}

simu.random <- function(p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list){
  idx <- sample(1:length(U.test.list),0.4*length(U.test.list))
  stats <- matrix(0,1,length(U.test.list)*2)
  stats[1,idx] <- 1
  U <- matrix(0,p,p)
  for (i in 1:length(idx)){
    U <- U + U.test.list[[idx[i]]]
  }
  generation <- random.generation(p,n,U)
  X <- generation$X
  Y <- generation$Y
  stats[1,(length(U.test.list)+1):(length(U.test.list)*2)] <- simu.single.2lr(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list)
  return(stats)
}

simu.random2 <- function(p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list){
  idx <- sample(1:length(U.test.list),0.4*length(U.test.list))
  stats <- matrix(0,1,length(U.test.list)*2)
  stats[1,idx] <- 1
  U <- matrix(0,p,p)
  for (i in 1:length(idx)){
    U <- U + U.test.list[[idx[i]]]
  }
  generation <- random.generation2(p,n,U)
  X <- generation$X
  Y <- generation$Y
  stats[1,(length(U.test.list)+1):(length(U.test.list)*2)] <- simu.single.2lr(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list)
  return(stats)
}

clusterExport(cl,c("simu.single.2lr","is.acyclic","random.generation","random.generation2","simu.random","simu.random2"))

## random graph 1
simu <- function(i,U.test.list){
  p <- 100;n <- 500;
  tau.list1 <- seq(0.1,0.2,0.01)
  gamma.list1 <- seq(0.05,0.5,0.05)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  tau.list3 <- seq(0.05,0.1,0.01)
  gamma.list3 <- seq(0.005,0.1,0.005)
  n.fold1 <- n.fold2 <- n.fold3 <- 5
  stats.random <- numeric(100)
  stats.random <- simu.random(p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list)
  return(stats.random)
}

set.seed(1); U.test.list <- test.generation(100,50)
len.test <- 100
clusterExport(cl,c("simu","len.test","U.test.list"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i,U.test.list),stat<-rep(-1,len.test));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=100)
write.csv(stats_random,file.path("./primary_results/part5/","stats_random1.csv"),row.names = FALSE)

## random graph 2
simu <- function(i,U.test.list){
  p <- 100;n <- 400;
  tau.list1 <- seq(0.1,0.2,0.01)
  gamma.list1 <- seq(0.05,0.5,0.05)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  tau.list3 <- seq(0.05,0.1,0.01)
  gamma.list3 <- seq(0.005,0.1,0.005)
  n.fold1 <- n.fold2 <- n.fold3 <- 5
  stats.random <- numeric(100)
  stats.random <- simu.random(p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list)
  return(stats.random)
}

set.seed(2); U.test.list <- test.generation(100,50)
len.test <- 100
clusterExport(cl,c("simu","len.test","U.test.list"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i,U.test.list),stat<-rep(-1,len.test));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=100)
write.csv(stats_random,file.path("./primary_results/part5/","stats_random2.csv"),row.names = FALSE)

## random graph 3
simu <- function(i,U.test.list){
  p <- 100;n <- 300;
  tau.list1 <- seq(0.1,0.2,0.01)
  gamma.list1 <- seq(0.05,0.5,0.05)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  tau.list3 <- seq(0.05,0.1,0.01)
  gamma.list3 <- seq(0.005,0.1,0.005)
  n.fold1 <- n.fold2 <- n.fold3 <- 5
  stats.random <- numeric(100)
  stats.random <- simu.random(p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list)
  return(stats.random)
}

set.seed(3); U.test.list <- test.generation(100,50)
len.test <- 100
clusterExport(cl,c("simu","len.test","U.test.list"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i,U.test.list),stat<-rep(-1,len.test));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=100)
write.csv(stats_random,file.path("./primary_results/part5/","stats_random3.csv"),row.names = FALSE)

## random graph 4
simu <- function(i,U.test.list){
  p <- 100;n <- 500;
  tau.list1 <- seq(0.1,0.2,0.01)
  gamma.list1 <- seq(0.05,0.5,0.05)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  tau.list3 <- seq(0.05,0.1,0.01)
  gamma.list3 <- seq(0.005,0.1,0.005)
  n.fold1 <- n.fold2 <- n.fold3 <- 5
  stats.random <- numeric(100)
  stats.random <- simu.random2(p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list)
  return(stats.random)
}

set.seed(4); U.test.list <- test.generation(100,50)
len.test <- 100
clusterExport(cl,c("simu","len.test","U.test.list"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i,U.test.list),stat<-rep(-1,len.test));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=100)
write.csv(stats_random,file.path("./primary_results/part5/","stats_random4.csv"),row.names = FALSE)

## random graph 5
simu <- function(i,U.test.list){
  p <- 100;n <- 400;
  tau.list1 <- seq(0.1,0.2,0.01)
  gamma.list1 <- seq(0.05,0.5,0.05)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  tau.list3 <- seq(0.05,0.1,0.01)
  gamma.list3 <- seq(0.005,0.1,0.005)
  n.fold1 <- n.fold2 <- n.fold3 <- 5
  stats.random <- numeric(100)
  stats.random <- simu.random2(p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list)
  return(stats.random)
}

set.seed(5); U.test.list <- test.generation(100,50)
len.test <- 100
clusterExport(cl,c("simu","len.test","U.test.list"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i,U.test.list),stat<-rep(-1,len.test));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=100)
write.csv(stats_random,file.path("./primary_results/part5/","stats_random5.csv"),row.names = FALSE)

## random graph 6
simu <- function(i,U.test.list){
  p <- 100;n <- 300;
  tau.list1 <- seq(0.1,0.2,0.01)
  gamma.list1 <- seq(0.05,0.5,0.05)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  tau.list3 <- seq(0.05,0.1,0.01)
  gamma.list3 <- seq(0.005,0.1,0.005)
  n.fold1 <- n.fold2 <- n.fold3 <- 5
  stats.random <- numeric(100)
  stats.random <- simu.random2(p,n,tau.list1,gamma.list1,tau.list2,gamma.list2,tau.list3,gamma.list3,n.fold1,n.fold2,n.fold3,U.test.list)
  return(stats.random)
}

set.seed(6); U.test.list <- test.generation(100,50)
len.test <- 100
clusterExport(cl,c("simu","len.test","U.test.list"))
stats.random <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i,U.test.list),stat<-rep(-1,len.test));return(stat)})
stats_random <- matrix(unlist(stats.random),byrow = TRUE,ncol=100)
write.csv(stats_random,file.path("./primary_results/part5/","stats_random6.csv"),row.names = FALSE)


stopCluster(cl)

