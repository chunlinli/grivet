library(lrpsadmm)
library(mvtnorm)
library(bnlearn)
library(pcalg)
library(grivet)

## This function implements LRpS+GES, for details, refer to https://github.com/benjaminfrot/lrpsadmm-examples/blob/master/LRpSGES.ipynb
causal_discovery_LRpS_GES <- function(y) {
  y <- scale(y)
  n <- nrow(y)
  p <- ncol(y)
  gammas <- c(0.05, 0.07, 0.1, 0.12, 0.15, 0.17, 0.2)
  ## Compute LRpS and select a graphical model using cross-validation
  xval_path <- lrpsadmm.cv(
    X = y,
    gammas = gammas,
    covariance.estimator = cor,
    n.folds = 5,
    verbose = FALSE,
    n.lambdas = 40,
    lambda.ratio = 1e-04,
    backend = "RcppEigen"
  )
  # Because the GES function of pcalg can only take a data matrix as input and not a covariance matrix, 
  # we simulate data with the exact same sample covariance matrix as the estimated one. 
  # This is a trick that lets use the GES function of pcalg.
  selected_s <- xval_path$best.fit$fit$S
  fake_data <- generate.data.for.GES(
    Sest = selected_s,
    n = n,
    p = p
  )
  # We now compute the GES path on the resulting data and compute the BIC at each step. 
  lrps_ges_output <- run.GES.and.select.with.BIC(
    obs.data = fake_data,
    nv = p,
    sim.data = data.frame(y)
  )
  
  list(out = lrps_ges_output, u = lrps_ges_output$best.essgraph)
}

## This function implements RFCI algorithm
causal_discovery_RFCI <- function(y) {
  suffStat <- list(C = cor(y), n = nrow(y))
  out <- rfci(
    suffStat = suffStat,
    p = ncol(y),
    skel.method = "stable",
    indepTest = gaussCItest,
    alpha = 0.001,
    verbose = FALSE
  )
  
  list(out = out, u = out@amat)
}

## This function estimates the structure of DAG using our prposed GrIVET method.

## Arguments:
## X: n*q data matrix for intervention variables.
## Y: n*p data matrix for primary variables.
## tau.list1 & gamma.list1: lists of tunning parameters for TLP regression in matrix V estimation.
## tau.list2 & gamma.list2; lists of tunning parameters for TLP regression in structure learning(penalized coefficient estimation).
## n.fold1 & n.fold2: number of folds for cross-validationg in tunning parameters selection, respectively for V estimation and structure learning.

## Returns:
## out: a p*p matrix representing the structure learning results, 1 for existing edge and 0 if non-existing.
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

## This function transmit a CPDAG into a DAG with a favorable fashion
PAG2best <- function(true_graph,estimate_graph){
  true_graph <- (true_graph != 0) * 1
  u1 <- pag2edge(estimate_graph)+1
  u2 <- 1*(u1==2)+1*(u1==1)*(true_graph==1)
  return(u2)
}

## This function evaluates the performance of structure learning.

## Arguments:
## true_graph: a p*p matrix representing the underlying structure of DAG, 1 if edges exist and 0 otherwise.
## estimate_graph: a p*p matrix representing the estimated structure of DAG, 1 if edges estimated to be existent and 0 otherwise.

## Returns:
## FDR: false discovery rate.
## SHD: structural hamming distance.
## JCI: jaccard index.
## TPR: true positive rate.
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

## The function `generate.data.for.GES` is from https://github.com/benjaminfrot/lrpsadmm-examples/blob/master/utils/generate_data_for_GES.R

generate.data.for.GES <- function(Sest, n, p) {
  S.est.lrps <- Sest
  Sig.est.lrps <- solve(S.est.lrps)
  fake.data <- rmvnorm(max(n,2*p), mean = rep(0, p), sigma = Sig.est.lrps)
  e <- eigen(Sig.est.lrps)
  sqrt.true.cov.mat <- e$vectors%*%sqrt(diag(e$values))
  samp.cov.mat <- cov(fake.data)
  e <- eigen(samp.cov.mat)
  sqrt.samp.cov.mat <- e$vectors%*%sqrt(diag(e$values))
  fake.data <- t(sqrt.true.cov.mat%*%solve(sqrt.samp.cov.mat,t(fake.data)))
  fake.data <- as.data.frame(fake.data)
  obs.data <- fake.data
}


## The function `run.GES.and.select.with.BIC` is from: https://github.com/benjaminfrot/lrpsadmm-examples/blob/master/utils/run_GES.R

run.GES.and.select.with.BIC <- function(obs.data, nv, sim.data) {
  
  rho <- 100000 # Compute the path, starting with this value of Rho
  rho.base <- 1.1 # The next value of Rho is Rho / 1.1
  path <- list()
  counter <- 1
  while(T) {
    score <- new("GaussL0penObsScore", obs.data, lambda = rho)
    start.time <- Sys.time()
    ges.fit <- ges(score)
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    
    u.amat <- as(ges.fit$essgraph, "matrix") * 1
    if(sum(u.amat) == 0) {
      rho <- rho / rho.base
      next()
    }
    
    # Compute the BIC
    bn <- empty.graph(colnames(obs.data))
    amat(bn) <- as(ges.fit$repr, "matrix") 
    bn <- bn.fit(bn, as.data.frame(obs.data))
    BIC <- BIC(bn, data = as.data.frame(obs.data))
    LogLik <- logLik(bn, data = as.data.frame(obs.data))
    NEdges <- sum(amat(bn)!=0)
    
    est.cpdag <- as(ges.fit$essgraph, "matrix") * 1
    # Compare to the true CPDAG
    #perf.metrics <- compute_metrics(true.dag = sim.data$true.obs.dag.amat, 
    #                                est.cpdag = est.cpdag)
    
    # Record these results
    path[[counter]] <- list()
    path[[counter]]$rho <- rho
    #path[[counter]]$metric <- perf.metrics
    path[[counter]]$NEdges <- NEdges
    path[[counter]]$BIC <- BIC
    path[[counter]]$LogLik <- LogLik
    path[[counter]]$fitting.time <- time.taken
    
    if (NEdges / choose(nv, 2) > 0.5) {
      break()
    }
    
    counter <- counter + 1
    rho <- rho / rho.base
  }
  
  # Get the value of lambda that gives the best BIC
  BIC <- unlist(sapply(path, function(a){get("BIC", a)}))
  idx <- which.max(BIC)
  rho.bic <- unlist(sapply(path, function(a) {get("rho", a)}))[idx]
  
  # Refit the model with this value of Rho
  score <- new("GaussL0penObsScore", obs.data, lambda = rho.bic)
  ges.fit <- ges(score)
  u.amat <- as(ges.fit$essgraph, "matrix") * 1
  
  list(path=path, best.essgraph=u.amat, best.fit=ges.fit)
}

## This function generate a hub graph with continuous intervention variables..

## Arguments: 
## p: the number of primary variables.
## n: the number of observations.

## Returns:
## X: a n*q matrix for intervention variables.
## Y: a n*p matrix for primary variables.
## U: a p*P matrix representing the caual effect matrix.
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

## This function generate a hub graph with discrete intervention variables..

## Arguments: 
## p: the number of primary variables.
## n: the number of observations.

## Returns:
## X: a n*q matrix for intervention variables.
## Y: a n*p matrix for primary variables.
## U: a p*P matrix representing the caual effect matrix.
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

## This function generate a random graph with continuous intervention variables..

## Arguments: 
## p: the number of primary variables.
## n: the number of observations.
## U.test: a p*p matrix, 1 representing tested edges, 0 otherwise.
## type: type = 0 controls the simulated graph so that H_0 holds, otherwise H_1 holds.

## Returns:
## X: a n*q matrix for intervention variables.
## Y: a n*p matrix for primary variables.
## U: a p*P matrix representing the caual effect matrix.
## W: a q*p matrix representing the interventional effect matrix.
## Sigma: a p*p matrix representing the covariance matrix of residuals.
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

## This function generate a random graph with discrete intervention variables..

## Arguments: 
## p: the number of primary variables.
## n: the number of observations.
## U.test: a p*p matrix, 1 representing tested edges, 0 otherwise.
## type: type = 0 controls the simulated graph so that H_0 holds, otherwise H_1 holds.

## Returns:
## X: a n*q matrix for intervention variables.
## Y: a n*p matrix for primary variables.
## U: a p*P matrix representing the caual effect matrix.
## W: a q*p matrix representing the interventional effect matrix.
## Sigma: a p*p matrix representing the covariance matrix of residuals.

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
