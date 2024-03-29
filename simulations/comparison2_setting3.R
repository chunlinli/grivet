library(grivet)

## This function estimates the coefficient without taking considerations of confounders.

## Arguments:
## X: a n*q matrix for intervention variables.
## Y: a n*p matrix for primary variables.
## Pi: the matrix representing the support of ancestral relationships, 1 for existing ancestral relation, 0 otherwise.
## Phi: the matrix representing the support of interventional relationships, 1 for existing interventional relation, 0 otherwise.
## tau.list & gamma.list: a list of tuning parameters to be selected for TLP regression.
## n.fold: the number of folds in selecting tuning parameters for TLP regression.

## Returns:
## U: the estimated causal effect matrix.
coef.direct <- function(X,Y,Pi,Phi,tau.list,gamma.list,n.fold){
  q <- ncol(X)
  p <- ncol(Y)
  U <- matrix(0,p,p)
  
  ind1 <- which(colMeans(Pi)!=0)
  if (length(ind1)==0){
    return(U)
  }
  
  for (i in 1:length(ind1)){
    child.idx <- ind1[i]
    Y.idx <- which(Pi[,child.idx]==1)
    X.idx <- which(Phi[,child.idx]==1)
    Z <- cbind(Y[,Y.idx],X[,X.idx])
    y <- Y[,child.idx]
    res <- cv.tlpreg0.aic(y,Z,tau.list=tau.list,gamma.list=gamma.list,n.fold=n.fold)
    tau <- res$tau
    gamma <- res$gamma
    b <- tlpreg0(y,Z,tau=tau,gamma=gamma)
    ind <- b!=0
    b.ols <- matrix(0,ncol(Z),1)
    if (sum(ind)==0){
      b.ols[,1] <- b
    }else{
      aic.memo <- numeric(sum(ind))
      ord <- order(abs(b),decreasing = TRUE)
      for (l in 1:sum(ind)){
        ind2 <- which(abs(b)>=abs(b[ord[l]]))
        mod <- lm(y~Z[,ind2,drop=FALSE])
        aic.memo[l] <- AIC(mod)
      }
      l <- which.min(aic.memo)
      ind <- which(abs(b)>=abs(b[ord[l]]))
      mod <- lm(y~Z[,ind,drop=FALSE])
      b.ols[ind,1] <- as.numeric(mod$coefficients)[-1]
    }
    
    U[Y.idx,child.idx] <- b.ols[1:length(Y.idx),1]
  }
  return(U)
}

## This function apply our GrIVET to estimate the coefficients.

## Arguments:
## X: a n*q matrix for intervention variables.
## Y: a n*p matrix for primary variables.
## Pi: the matrix representing the support of ancestral relationships, 1 for existing ancestral relation, 0 otherwise.
## Phi: the matrix representing the support of interventional relationships, 1 for existing interventional relation, 0 otherwise.
## Piv: the matrix representing the candidate IV set, 1 for existing candidate IV relationship and 0 otherwise.
## tau.list & gamma.list: a list of tuning parameters to be selected for TLP regression.
## n.fold: the number of folds in selecting tuning parameters for TLP regression.

## Returns:
## U: the estimated causal effect matrix.
coef.proposed <- function(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold){
  out <- cv.intdag.coe(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold)
  U <- out$U
  return(U)
}

## This function evaluate the performance of coefficient estimation.

## Arguments:
## U: the true causal effect matrix.
## U_es: the estimated causal effect matrix.
## Pi: the matrix representing true support of ancestral relationships, 1 for existing ancestral relation, 0 otherwise.

## Returns:
## max_abs_diff: maximum absolute difference between U and U_es.
## mean_abs_diff: mean absolute difference between U and U_es.
## mean_sq_diff: mean square difference between U and U_es.
metrics <- function(U,U_es,Pi){
  max_abs_diff <- max(abs(U-U_es))
  mean_abs_diff <- sum(abs(U-U_es))/max(1,sum(Pi!=0))
  mean_sq_diff <- sum((U-U_es)^2)/max(1,sum(Pi!=0))
  return(c(max_abs_diff, mean_abs_diff,mean_sq_diff))
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
  vari <- diag(runif(p,min=0.8,max=1.2))
  sigma <- vari%*%diag(p)%*%vari
  e <- mvrnorm(n,rep(0,p),sigma)
  r <- 2
  V <- matrix(0,r,p)
  V[1,1:(p/2)] <- 1; V[2,(p/2 + 1):p] <- 1
  V <- V*matrix(runif(r*p,min=0.8,max=1.2)*sample(c(-1,1),r*p,replace = TRUE),r,p)
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
  vari <- diag(runif(p,min=0.8,max=1.2))
  sigma <- vari%*%diag(p)%*%vari
  e <- mvrnorm(n,rep(0,p),sigma)
  r <- 2
  V <- matrix(0,r,p)
  V[1,1:(p/2)] <- 1; V[2,(p/2 + 1):p] <- 1
  V <- V*matrix(runif(r*p,min=0.8,max=1.2)*sample(c(-1,1),r*p,replace = TRUE),r,p)
  Z <- matrix(rnorm(r*n),ncol=r)
  X <- matrix(sample(c(-1,1),q*n,replace = TRUE),ncol=q); Y <- matrix(0,n,p); Y <- (X%*%W+e+Z%*%V)%*%solve(diag(1,p)-U);
  Sigma <- sigma+t(V)%*%V
  return(list(X=X,Y=Y,U=U,W=W,Sigma=Sigma))
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
  vari <- diag(runif(p,min=0.8,max=1.2))
  sigma <- vari%*%diag(p)%*%vari
  e <- mvrnorm(n,rep(0,p),sigma)
  r <- 2
  V <- matrix(0,r,p)
  V[1,1:((p+1)/2)] <- 1; V[2,((p+3)/2):p] <- 1;
  V <- V*matrix(runif(r*p,min=0.8,max=1.2)*sample(c(-1,1),r*p,replace = TRUE),r,p)
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
  vari <- diag(runif(p,min=0.8,max=1.2))
  sigma <- vari%*%diag(p)%*%vari
  e <- mvrnorm(n,rep(0,p),sigma)
  r <- 2
  V <- matrix(0,r,p)
  V[1,1:((p+1)/2)] <- 1; V[2,((p+3)/2):p] <- 1;
  V <- V*matrix(runif(r*p,min=0.8,max=1.2)*sample(c(-1,1),r*p,replace = TRUE),r,p)
  Z <- matrix(rnorm(r*n),ncol=r)
  X <- matrix(sample(c(-1,1),q*n,replace = TRUE),ncol=q); Y <- matrix(0,n,p); Y <- (X%*%W+e+Z%*%V)%*%solve(diag(1,p)-U);
  Sigma <- sigma+t(V)%*%V
  return(list(X=X,Y=Y,U=U,W=W,Sigma=Sigma))
}