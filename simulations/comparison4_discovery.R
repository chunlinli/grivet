library(doParallel)
n_cores <- detectCores() - 1
registerDoParallel(cores=n_cores)  
cl <- makeCluster(n_cores, type="FORK")

clusterEvalQ(cl,{library(grivet)
  library(mvtnorm)
  library(MASS)})

causal_discovery_proposed <- function(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2){
  result.stage1_1 <- cv.intdag.pmle.diff.aic(X,Y,tau.list1,gamma.list1,n.fold1)
  V.es <- result.stage1_1$V
  result.stage1_2 <- topological_order(V.es)
  Pi_es <- result.stage1_2$an_mat
  Phi_es <- result.stage1_2$in_mat
  Piv_es <- result.stage1_2$iv_mat
  result.stage2 <- cv.intdag.coe.penalized(X,Y,Pi_es,Phi_es,Piv_es,tau.list2,gamma.list2,n.fold2)
  out <- 1*(result.stage2$U!=0)
  return(out)
}

metrics <- function(true_graph, estimate_graph) {
  
  true_edges <- (true_graph != 0) * 1
  skeleton <- (t(true_graph) != 0) * 1 + (true_graph != 0) * 1
  undirected_edges <- (estimate_graph*t(estimate_graph))
  directed_edges <- estimate_graph - undirected_edges
  
  # true positive
  true_positive_directed <- directed_edges * true_edges
  true_positive_undirected <- undirected_edges * true_edges
  true_positive <- true_positive_directed + true_positive_undirected
  
  # false positive
  false_positive_directed <- directed_edges * (true_graph != 1)
  false_positive_undirected <- undirected_edges * (skeleton != 1)
  false_positive <- false_positive_directed + false_positive_undirected
  
  # missing
  missing <- true_edges - true_positive
  
  estimated_size <- sum(directed_edges) + sum(undirected_edges)/2 
  
  FDR <- (sum(false_positive_directed)+sum(false_positive_undirected)/2) / max(1, estimated_size)
  SHD <- sum(false_positive_directed) + sum(missing) + sum(false_positive_undirected)/2
  JCI <- sum(true_positive) / max(1,sum(true_positive)+SHD)
  TPR <- sum(true_positive) / max(1, sum(true_edges))
  return(c(FDR,JCI,TPR,SHD))
}

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

clusterExport(cl,c("causal_discovery_proposed","metrics","hub.generation","hub.generation2"))

## hub graph 1
simu <- function(i){
  p <- 10; inv <- 2; n <- 1000;
  tau.list1 <- seq(0.2,0.3,0.01)
  gamma.list1 <- seq(0.1,1,0.1)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  n.fold1 <- n.fold2 <- 5
  res <- numeric(24)
  
  for (q in 12:17){
    generation <- hub.generation(p,q,inv,n)
    X <- generation$X
    Y <- generation$Y
    U <- generation$U
    U <- 1*(U!=0)
    u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
    res[((q-12)*4+1):((q-11)*4)] <- metrics(U,u3)
  }
  
  return(res)
}

len.stat <- 24
clusterExport(cl,c("simu","len.stat"))
stats.hub <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_hub <- matrix(unlist(stats.hub),byrow = TRUE,ncol=len.stat)
write.csv(stats_hub,file.path("./primary_results/part4/discovery/","stats_hub1.csv"),row.names = FALSE)

## hub graph 2
simu <- function(i){
  p <- 10; inv <- 2; n <- 500;
  tau.list1 <- seq(0.2,0.3,0.01)
  gamma.list1 <- seq(0.1,1,0.1)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  n.fold1 <- n.fold2 <- 5
  res <- numeric(24)
  
  for (q in 12:17){
    generation <- hub.generation(p,q,inv,n)
    X <- generation$X
    Y <- generation$Y
    U <- generation$U
    U <- 1*(U!=0)
    u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
    res[((q-12)*4+1):((q-11)*4)] <- metrics(U,u3)
  }
  
  return(res)
}

len.stat <- 24
clusterExport(cl,c("simu","len.stat"))
stats.hub <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_hub <- matrix(unlist(stats.hub),byrow = TRUE,ncol=len.stat)
write.csv(stats_hub,file.path("./primary_results/part4/discovery/","stats_hub2.csv"),row.names = FALSE)

## hub graph 3
simu <- function(i){
  p <- 10; inv <- 2; n <- 200;
  tau.list1 <- seq(0.2,0.3,0.01)
  gamma.list1 <- seq(0.1,1,0.1)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  n.fold1 <- n.fold2 <- 5
  res <- numeric(24)
  
  for (q in 12:17){
    generation <- hub.generation(p,q,inv,n)
    X <- generation$X
    Y <- generation$Y
    U <- generation$U
    U <- 1*(U!=0)
    u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
    res[((q-12)*4+1):((q-11)*4)] <- metrics(U,u3)
  }
  
  return(res)
}

len.stat <- 24
clusterExport(cl,c("simu","len.stat"))
stats.hub <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_hub <- matrix(unlist(stats.hub),byrow = TRUE,ncol=len.stat)
write.csv(stats_hub,file.path("./primary_results/part4/discovery/","stats_hub3.csv"),row.names = FALSE)

## hub graph 4
simu <- function(i){
  p <- 10; inv <- 2; n <- 1000;
  tau.list1 <- seq(0.2,0.3,0.01)
  gamma.list1 <- seq(0.1,1,0.1)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  n.fold1 <- n.fold2 <- 5
  res <- numeric(24)
  
  for (q in 12:17){
    generation <- hub.generation2(p,q,inv,n)
    X <- generation$X
    Y <- generation$Y
    U <- generation$U
    U <- 1*(U!=0)
    u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
    res[((q-12)*4+1):((q-11)*4)] <- metrics(U,u3)
  }
  
  return(res)
}

len.stat <- 24
clusterExport(cl,c("simu","len.stat"))
stats.hub <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_hub <- matrix(unlist(stats.hub),byrow = TRUE,ncol=len.stat)
write.csv(stats_hub,file.path("./primary_results/part4/discovery/","stats_hub4.csv"),row.names = FALSE)

## hub graph 5
simu <- function(i){
  p <- 10; inv <- 2; n <- 500;
  tau.list1 <- seq(0.2,0.3,0.01)
  gamma.list1 <- seq(0.1,1,0.1)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  n.fold1 <- n.fold2 <- 5
  res <- numeric(24)
  
  for (q in 12:17){
    generation <- hub.generation2(p,q,inv,n)
    X <- generation$X
    Y <- generation$Y
    U <- generation$U
    U <- 1*(U!=0)
    u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
    res[((q-12)*4+1):((q-11)*4)] <- metrics(U,u3)
  }
  
  return(res)
}

len.stat <- 24
clusterExport(cl,c("simu","len.stat"))
stats.hub <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_hub <- matrix(unlist(stats.hub),byrow = TRUE,ncol=len.stat)
write.csv(stats_hub,file.path("./primary_results/part4/discovery/","stats_hub5.csv"),row.names = FALSE)

## hub graph 6
simu <- function(i){
  p <- 10; inv <- 2; n <- 200;
  tau.list1 <- seq(0.2,0.3,0.01)
  gamma.list1 <- seq(0.1,1,0.1)
  tau.list2 <- seq(0.1,0.2,0.01)
  gamma.list2 <- seq(0.05,0.5,0.05)
  n.fold1 <- n.fold2 <- 5
  res <- numeric(24)
  
  for (q in 12:17){
    generation <- hub.generation2(p,q,inv,n)
    X <- generation$X
    Y <- generation$Y
    U <- generation$U
    U <- 1*(U!=0)
    u3 <- causal_discovery_proposed(X,Y,tau.list1,gamma.list1,tau.list2,gamma.list2,n.fold1,n.fold2)
    res[((q-12)*4+1):((q-11)*4)] <- metrics(U,u3)
  }
  
  return(res)
}

len.stat <- 24
clusterExport(cl,c("simu","len.stat"))
stats.hub <- parLapply(cl,1:1000,function(i){set.seed(i);try(stat<-simu(i),stat<-rep(-1,len.stat));return(stat)})
stats_hub <- matrix(unlist(stats.hub),byrow = TRUE,ncol=len.stat)
write.csv(stats_hub,file.path("./primary_results/part4/discovery/","stats_hub6.csv"),row.names = FALSE)

stopCluster(cl)
