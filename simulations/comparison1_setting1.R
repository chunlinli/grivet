library(mvtnorm)
library(grivet)

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

hub.generation <- function(p,n){
  if ((p-1)%%10!=0)
    stop("wrong p is given")
  q <- 2*p+(p-1)/2
  U <- matrix(0,p,p);U[1,2:p] <- 1;
  W1 <- W2 <- diag(1,p,p); W3 <- matrix(0,q-2*p,p); 
  for (i in 1:(q-2*p)){
    W3[i,2*i]<-1
    W3[i,(2*i+1)]<-1
  }
  W <- rbind(W1,W2,W3)
  U <- U*matrix(rep(1,p*p)*sample(c(-1,1),p*p,replace = TRUE),p,p)
  W <- W*matrix(rep(1,q*p)*sample(c(-1,1),q*p,replace = TRUE),q,p)
  vari <- diag(runif(p,min=0.4,max=0.6))
  sigma <- vari%*%diag(p)%*%vari
  e <- mvrnorm(n,rep(0,p),sigma)
  r <- (p-1)/10
  V <- matrix(0,r,p)
  for (i in 1:r){
    V[i,(10*(i-1)+2):(10*i+1)] <- 1
  }
  V[1,1] <- 1
  V <- V*matrix(runif(r*p,min=0.4,max=0.6)*sample(c(-1,1),r*p,replace = TRUE),r,p)
  Z <- matrix(rnorm(r*n),ncol=r)
  X <- matrix(rnorm(q*n),ncol=q); Y <- matrix(0,n,p); Y <- (X%*%W+e+Z%*%V)%*%solve(diag(1,p)-U);
  Sigma <- sigma+t(V)%*%V
  return(list(X=X,Y=Y,U=U,W=W,Sigma=Sigma))
}


hub.generation2 <- function(p,n){
  if ((p-1)%%10!=0)
    stop("wrong p is given")
  q <- 2*p+(p-1)/2
  U <- matrix(0,p,p);U[1,2:p] <- 1;
  W1 <- W2 <- diag(1,p,p); W3 <- matrix(0,q-2*p,p); 
  for (i in 1:(q-2*p)){
    W3[i,2*i]<-1
    W3[i,(2*i+1)]<-1
  }
  W <- rbind(W1,W2,W3)
  U <- U*matrix(rep(1,p*p)*sample(c(-1,1),p*p,replace = TRUE),p,p)
  W <- W*matrix(rep(1,q*p)*sample(c(-1,1),q*p,replace = TRUE),q,p)
  vari <- diag(runif(p,min=0.4,max=0.6))
  sigma <- vari%*%diag(p)%*%vari
  e <- mvrnorm(n,rep(0,p),sigma)
  r <- (p-1)/10
  V <- matrix(0,r,p)
  for (i in 1:r){
    V[i,(10*(i-1)+2):(10*i+1)] <- 1
  }
  V[1,1] <- 1
  V <- V*matrix(runif(r*p,min=0.4,max=0.6)*sample(c(-1,1),r*p,replace = TRUE),r,p)
  Z <- matrix(rnorm(r*n),ncol=r)
  X <- matrix(sample(c(-1,1),q*n,replace = TRUE),ncol=q); Y <- matrix(0,n,p); Y <- (X%*%W+e+Z%*%V)%*%solve(diag(1,p)-U);
  Sigma <- sigma+t(V)%*%V
  return(list(X=X,Y=Y,U=U,W=W,Sigma=Sigma))
}

random.generation <- function(p,n,U.test,type){
  if (p%%10!=0)
    stop("wrong p is given")
  q <- 2.5*p
  U <- matrix(sample(x=c(1,0),size=p*p,replace=TRUE,prob=c(1/(10*p),1-1/(10*p))),p,p)
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

random.generation2 <- function(p,n,U.test,type){
  if (p%%10!=0)
    stop("wrong p is given")
  q <- 2.5*p
  U <- matrix(sample(x=c(1,0),size=p*p,replace=TRUE,prob=c(1/(10*p),1-1/(10*p))),p,p)
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
