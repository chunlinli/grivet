## This function computes U and W that maximize the likelihood with a closed-form formula given the nonzero indices.
## The computational formula can be found in Appendix.
## X: A n*q matrix representing interventions.
## Y: A n*p matrix representing primary variables.
## U0: Matrix of Indices for non-zero elements of U.
## W0: Matrix of Indices for non-zero elements of W.
## S: Precision matrix.

intdag.mle <- function(X,Y,U0,W0,S){
  
  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(Y)
  V0 <- rbind(U0,W0)
  Z <- cbind(Y,X)
  Z2 <- crossprod(Y,Z)
  
  A1 <- crossprod(Z)
  A2 <- kronecker(S,A1)
  idx <- (as.vector(V0) == 1)
  A <- A2[idx,idx]
  
  B1 <- as.vector(crossprod(Z2,S))
  B <- B1[idx]
  
  try(res <- solve(A,B),{res <- solve(A+1e-6*diag(sum(idx)),B);cat("GG")})
  
  res1 <- numeric(p*(p+q))
  res1[idx] <- res
  res.mat <- matrix(res1,ncol = p,byrow = FALSE)
  
  U <- res.mat[1:p,]
  W <- res.mat[-(1:p),]
  return(list(U=U,W=W))
}

## This function computes U and W that maximize the likelihood with BCD algorithm given the nonzero indices.
## The computational algorithm can be found in Appendix.
## X: A n*q matrix representing interventions.
## Y: A n*p matrix representing primary variables.
## U0: Matrix of Indices for non-zero elements of U.
## W0: Matrix of Indices for non-zero elements of W.
## U.init: Initial value for U
## W.init: Initial value for W
## S: Precision matrix.
## max.it: Maximum number of iterations
## tol: Tolerance of stopping criterion. 

intdag.mle.bcd <- function(X,Y,U0,W0,U.init,W.init,S,max.it,tol){
  
  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(Y)
  v <- p+q
  Z <- cbind(Y,X)
  
  V0 <- rbind(U0,W0)
  V.init <- rbind(U.init,W.init)
  idx <- (as.vector(V0)==1)
  ind <- which(idx)
  num.nonzero <- sum(idx)
  
  loc.mat <- matrix(0,num.nonzero,2)
  for(i in 1:num.nonzero){
    ind.curr <- ind[i]
    loc.mat[i,2] <- ceiling(ind.curr/v)
    loc.mat[i,1] <- ind.curr - (loc.mat[i,2]-1)*v
  }
  
  Z.init <- Y - Z%*%V.init
  Z.prod <- crossprod(Z.init)
  T.pre <- T.curr <- sum(diag(Z.prod%*%S))
  tol.iter <- 1
  iter <- 1
  
  V.prev <- V.curr <- V.init
  
  
  while(tol.iter > tol&iter<max.it){
    for(i in 1:num.nonzero){
      Z.curr <- Y - Z%*%V.curr
      Z.prod <- crossprod(Z.curr)
      omicron.c <- sum(diag(Z.prod%*%S))
      y <- Z[,loc.mat[i,1],drop=FALSE]
      omicron.a <- S[loc.mat[i,2],loc.mat[i,2]]*crossprod(y)
      z <- crossprod(Z.curr,y)
      omicron.b <- -2*crossprod(z,S[,loc.mat[i,2]])
      if(omicron.a>0){
        omicron.increa <- -omicron.b/(2*omicron.a)
      }else{
        omicron.increa <- 0
      }
      V.curr[loc.mat[i,1],loc.mat[i,2]] <- V.curr[loc.mat[i,1],loc.mat[i,2]]+omicron.increa
    }
    
    Z.curr <- Y - Z%*%V.curr
    Z.prod <- crossprod(Z.curr)
    T.curr <- sum(diag(Z.prod%*%S))
    tol.iter <- T.pre - T.curr
    if(tol.iter < 0){
      V.curr <- V.prev
      break
    }
    V.prev <- V.curr
    T.pre <- T.curr
    
    iter <- iter+1
  }
  
  if(iter==max.it){
    warning("Coordinate Descent Algorithm for statistic computing does not converge!")
  }
  U <- V.curr[1:p,]
  W <- V.curr[-(1:p),]
  return(list(U=U,W=W))
}


#' @title Statistic Computation
#' @description This function computes 2Lr on the whole graph.
#' @param X A n*q matrix representing interventions.
#' @param Y A n*p matrix representing primary variables.
#' @param U0 Matrix of Indices for non-zero elements of U under the null hypothesis.
#' @param W0 Matrix of Indices for non-zero elements of W under the null hypothesis.
#' @param U1 Matrix of Indices for non-zero elements of U under the alternative hypothesis.
#' @param W1 Matrix of Indices for non-zero elements of W under the alternative hypothesis.
#' @param S Precision matrix.
#' @return A list consits of estimates of U,W under the null and alternative hypothesis and the statistic.
#' @references Discovery and inference of a causal network with hidden confounding.
#' @examples
#' \dontrun{
#' to be added.
#' }
#' @export
intdag.2lr <- function(X,Y,U0,W0,U1,W1,S){
  m0 <- intdag.mle(X,Y,U0,W0,S)
  U0.mle <- m0$U
  W0.mle <- m0$W
  Z0 <- Y - Y%*%U0.mle - X%*%W0.mle
  Z0 <- crossprod(Z0)
  T0 <- sum(diag(Z0%*%S))
  
  max.it <- 1000
  tol <- 0.01
  m1 <- intdag.mle.bcd(X,Y,U0,W0,U0.mle,W0.mle,S,max.it,tol)
  U1.mle <- m1$U
  W1.mle <- m1$W
  Z1 <- Y - Y%*%U1.mle - X%*%W1.mle
  Z1 <- crossprod(Z1)
  T1 <- sum(diag(Z1%*%S))
  
  m2 <- intdag.mle.bcd(X,Y,U1,W1,U1.mle,W1.mle,S,max.it,tol)
  U2.mle <- m2$U
  W2.mle <- m2$W
  Z2 <- Y - Y%*%U2.mle - X%*%W2.mle
  Z2 <- crossprod(Z2)
  T2 <- sum(diag(Z2%*%S))
  
  m3 <- intdag.mle(X,Y,U1,W1,S)
  U3.mle <- m3$U
  W3.mle <- m3$W
  Z3 <- Y - Y%*%U3.mle - X%*%W3.mle
  Z3 <- crossprod(Z3)
  T3 <- sum(diag(Z3%*%S))
  
  if(T3 < T2){
    m2 <- intdag.mle.bcd(X,Y,U1,W1,U3.mle,W3.mle,S,max.it,tol)
    U2.mle <- m2$U
    W2.mle <- m2$W
    Z2 <- Y - Y%*%U2.mle - X%*%W2.mle
    Z2 <- crossprod(Z2)
    T2 <- sum(diag(Z2%*%S))
  }
  
  t <- T1-T2
  
  return(list(U0.mle=U1.mle,W0.mle=W1.mle,U1.mle=U2.mle,W1.mle=W2.mle,statistic=t))
}


#' @title Statistic Computation
#' @description This function computes the local graph with nodes involved in the tests and their ancestors.
#' @param X A n*q matrix representing interventions.
#' @param Y A n*p matrix representing primary variables.
#' @param U.diff Edges to be tested in the hypothesis testing.
#' @param U.hat U or estimates of U.
#' @param W.hat W or estimates of W.
#' @param Pi Ancestral relationships.
#' @param Phi Interventional relationships.
#' @param Sigma Covariance matrix.
#' @return A list consits of related interventions X, related primary variable Y, and related estimates of U,W,Pi,Phi and Sigma.
#' @references Discovery and inference of a causal network with hidden confounding.
#' @examples
#' \dontrun{
#' to be added.
#' }
#' @export
intdag.localization.ver2 <- function(X,Y,U.diff,U.hat,W.hat,Pi,Phi,Sigma){
  Y.idx <- union(which(colSums(U.diff)!=0),which(rowSums(U.diff)!=0))
  Y.idx <- union(which(rowSums(Pi[,Y.idx])!=0),Y.idx)
  Y.idx <- Y.idx[order(Y.idx)]
  Pi.new <- Pi[Y.idx,Y.idx]
  U.test.new <- U.diff[Y.idx,Y.idx]
  Sigma.new <- Sigma[Y.idx,Y.idx]
  Y.new <- Y[,Y.idx,drop=FALSE]
  U.hat.new <- U.hat[Y.idx,Y.idx]
  Phi.new <- Phi[,Y.idx,drop=FALSE]
  X.idx <- which(rowSums(Phi.new)!=0)
  X.idx <- X.idx[order(X.idx)]
  Phi.new <- Phi.new[X.idx,,drop=FALSE]
  W.hat.new <- W.hat[X.idx,Y.idx,drop=FALSE]
  X.new <- X[,X.idx,drop=FALSE]
  return(list(X.new=X.new,Y.new=Y.new,U.test.new=U.test.new,U.hat.new=U.hat.new,W.hat.new=W.hat.new,Pi.new=Pi.new,Phi.new=Phi.new,Sigma.new=Sigma.new))  
}
