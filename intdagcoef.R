source("tlpreg.r",chdir = TRUE)
## This function implements coefficient estimation for stage 2. 
cv.intdag.coe <- function(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold){
  
  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(Y)
  
  Piv <- t(Piv)
  
  U <- matrix(0,p,p)
  W <- matrix(0,q,p)
  taus <- matrix(0,p,p)
  gammas <- matrix(0,p,p)
  
  ## estimate topological depth
  U.depth <- topo.depth(Pi)
  ## estimate computation orders
  U.order1 <- es.order(Pi,U.depth)
  U.order2 <- es.order.sub(U.order1)
  ## estimate relationship
  relations.list <- relations.finder(Pi,Phi)
  
  ## estimate U and W
  if(sum(U.order1 > 0)){
    U.order <- es.order.sub(U.order1)
    es.Y <- matrix(rep(0,n*p),nrow=n) ## es.Y represents estimates of primary variables
    
    ## Fill es.Y for primary variables with topological depth 0
    idx.0 <- which(U.depth==0)
    for (i in 1:length(idx.0)){
      idx.iter <- idx.0[i]
      idx.x <- which(Phi[,idx.iter]!=0)
      W[idx.x,idx.iter] <- as.numeric(solve(crossprod(X[,idx.x]),crossprod(X[,idx.x],Y[,idx.iter])))
      es.Y[,idx.iter] <- X[,idx.x,drop=FALSE]%*%W[idx.x,idx.iter]
    }
    
    for (layer.iter in 1:max(U.order1)){
      children.unique <- unique(U.order2[U.order2[,1]==layer.iter,3])
      for (child.iter in 1:length(children.unique)){
        child.idx <- children.unique[child.iter]
        parents.idx <- as.numeric(U.order2[U.order2[,1]==layer.iter & U.order2[,3]==child.idx,2])
        preparents.idx <- as.numeric(U.order2[U.order2[,1]<layer.iter & U.order2[,3]==child.idx,2])
        working.X <- es.Y[,parents.idx,drop=FALSE]
        if (length(preparents.idx) >0){
          working.Y <- Y[,child.idx] - Y[,preparents.idx,drop=FALSE]%*%U[preparents.idx,child.idx] 
        }else{
          working.Y <- Y[,child.idx,drop=FALSE] 
        }
        working.Z.idx1 <- which(colSums(relations.list[[2]][parents.idx,,drop=FALSE])>0)
        working.Z.idx2 <- which(colSums(relations.list[[4]][parents.idx,,drop=FALSE])+relations.list[[3]][child.idx,] >0)
        working.Z.idx3 <- which(colSums(Piv[parents.idx,,drop=FALSE])>0)
        working.Z.idx5 <- setdiff(working.Z.idx2,working.Z.idx3)
        working.Z.idx4 <- setdiff(working.Z.idx2,working.Z.idx5)
        if (length(working.Z.idx1)>0){
          working.Z1 <- cbind(Y[,working.Z.idx1,drop=FALSE],X[,working.Z.idx5,drop=FALSE])
        }else{
          working.Z1 <- X[,working.Z.idx5,drop=FALSE]
        }
        working.Z2 <- X[,working.Z.idx4,drop=FALSE]
        coef.iter.res <- coef.single(working.Y,working.X,working.Z1,working.Z2,tau.list,gamma.list,n.fold)
        coef.iter <- coef.iter.res$b
        taus[parents.idx,child.idx] <- coef.iter.res$tau
        gammas[parents.idx,child.idx] <- coef.iter.res$gamma
        U[parents.idx,child.idx] <- coef.iter
        
        ## update es.Y if necessary
        if(max(U.order1[,child.idx])<=layer.iter){
          ## updat W
          es.Y.part1 <- Y[,which(Pi[,child.idx]==1),drop=FALSE]%*%U[which(Pi[,child.idx]==1),child.idx]
          Y.resi <- Y[,child.idx] - es.Y.part1
          idx.x <- which(Phi[,child.idx]!=0)
          W[idx.x,child.idx] <- as.numeric(solve(crossprod(X[,idx.x]),crossprod(X[,idx.x],Y.resi)))
          ## update es.Y
          es.Y[,child.idx] <- es.Y.part1 + X[,idx.x,drop=FALSE]%*%W[idx.x,child.idx]
        }
      }
    }
  }else{
    ## In this situation, no links exist between primary variables, so we need only to estimate W
    for (j in 1:p){
      idx.x <- which(Phi[,j]!=0)
      W[idx.x,j]<- as.numeric(solve(crossprod(X[,idx.x]),crossprod(X[,idx.x],Y[,j])))
    }
    es.Y <- X%*%W
  }
  
  ## estimate the covariance matrix
  dfs <- colSums(Pi) + colSums(Phi)
  E <- Y - es.Y
  E.modified <- apply(E,2,function(e) e-mean(e))
  E.modified <- t(t(E.modified)/sqrt(n-sqrt(dfs)))
  Sigma <- crossprod(E.modified)
  
  return(list(U=U,W=W,Sigma=Sigma,taus=taus,gammas=gammas))
}

cv.intdag.coe.penalized <- function(X,Y,Pi,Phi,Piv,tau.list,gamma.list,n.fold){
  
  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(Y)
  
  Piv <- t(Piv)
  
  U <- matrix(0,p,p)
  W <- matrix(0,q,p)
  taus <- matrix(0,p,p)
  gammas <- matrix(0,p,p)
  
  ## estimate topological depth
  U.depth <- topo.depth(Pi)
  ## estimate computation orders
  U.order1 <- es.order(Pi,U.depth)
  U.order2 <- es.order.sub(U.order1)
  ## estimate relationship
  relations.list <- relations.finder(Pi,Phi)
  
  ## estimate U and W
  if(sum(U.order1 > 0)){
    U.order <- es.order.sub(U.order1)
    es.Y <- matrix(rep(0,n*p),nrow=n) ## es.Y represents estimates of primary variables
    
    ## Fill es.Y for primary variables with topological depth 0
    idx.0 <- which(U.depth==0)
    for (i in 1:length(idx.0)){
      idx.iter <- idx.0[i]
      idx.x <- which(Phi[,idx.iter]!=0)
      W[idx.x,idx.iter] <- as.numeric(solve(crossprod(X[,idx.x]),crossprod(X[,idx.x],Y[,idx.iter])))
      es.Y[,idx.iter] <- X[,idx.x,drop=FALSE]%*%W[idx.x,idx.iter]
    }
    
    for (layer.iter in 1:max(U.order1)){
      children.unique <- unique(U.order2[U.order2[,1]==layer.iter,3])
      for (child.iter in 1:length(children.unique)){
        child.idx <- children.unique[child.iter]
        parents.idx <- as.numeric(U.order2[U.order2[,1]==layer.iter & U.order2[,3]==child.idx,2])
        preparents.idx <- as.numeric(U.order2[U.order2[,1]<layer.iter & U.order2[,3]==child.idx,2])
        working.X <- es.Y[,parents.idx,drop=FALSE]
        if (length(preparents.idx) >0){
          working.Y <- Y[,child.idx] - Y[,preparents.idx,drop=FALSE]%*%U[preparents.idx,child.idx] 
        }else{
          working.Y <- Y[,child.idx,drop=FALSE] 
        }
        working.Z.idx1 <- which(colSums(relations.list[[2]][parents.idx,,drop=FALSE])>0)
        working.Z.idx2 <- which(colSums(relations.list[[4]][parents.idx,,drop=FALSE])+relations.list[[3]][child.idx,] >0)
        working.Z.idx3 <- which(colSums(Piv[parents.idx,,drop=FALSE])>0)
        working.Z.idx5 <- setdiff(working.Z.idx2,working.Z.idx3)
        working.Z.idx4 <- setdiff(working.Z.idx2,working.Z.idx5)
        if (length(working.Z.idx1)>0){
          working.Z1 <- cbind(Y[,working.Z.idx1,drop=FALSE],X[,working.Z.idx5,drop=FALSE])
        }else{
          working.Z1 <- X[,working.Z.idx5,drop=FALSE]
        }
        working.Z2 <- X[,working.Z.idx4,drop=FALSE]
        coef.iter.res <- coef.single.penalized(working.Y,working.X,working.Z1,working.Z2,tau.list,gamma.list,n.fold)
        coef.iter <- coef.iter.res$b
        taus[parents.idx,child.idx] <- coef.iter.res$tau
        gammas[parents.idx,child.idx] <- coef.iter.res$gamma
        U[parents.idx,child.idx] <- coef.iter
        
        ## update es.Y if necessary
        if(max(U.order1[,child.idx])<=layer.iter){
          ## updat W
          es.Y.part1 <- Y[,which(Pi[,child.idx]==1),drop=FALSE]%*%U[which(Pi[,child.idx]==1),child.idx]
          Y.resi <- Y[,child.idx] - es.Y.part1
          idx.x <- which(Phi[,child.idx]!=0)
          W[idx.x,child.idx] <- as.numeric(solve(crossprod(X[,idx.x]),crossprod(X[,idx.x],Y.resi)))
          ## update es.Y
          es.Y[,child.idx] <- es.Y.part1 + X[,idx.x,drop=FALSE]%*%W[idx.x,child.idx]
        }
      }
    }
  }else{
    ## In this situation, no links exist between primary variables, so we need only to estimate W
    for (j in 1:p){
      idx.x <- which(Phi[,j]!=0)
      W[idx.x,j]<- as.numeric(solve(crossprod(X[,idx.x]),crossprod(X[,idx.x],Y[,j])))
    }
    es.Y <- X%*%W
  }
  
  ## estimate the covariance matrix
  dfs <- colSums(Pi) + colSums(Phi)
  E <- Y - es.Y
  E.modified <- apply(E,2,function(e) e-mean(e))
  E.modified <- t(t(E.modified)/sqrt(n-sqrt(dfs)))
  Sigma <- crossprod(E.modified)
  
  return(list(U=U,W=W,Sigma=Sigma,taus=taus,gammas=gammas))
}

## This function returns estimate of a single element in matrix U
coef.single <- function(Y,X,Z1,Z2,tau.list,gamma.list,n.fold){
  X.cbind <- cbind(X,Z1,Z2)
  p <- ncol(X.cbind)
  q1 <- ncol(X)
  q2 <- ncol(Z1)
  pen.factors <- c(rep(0,q1),rep(1,p-q1))
  res <- cv.tlpreg0.aic(y=Y,X=X.cbind,tau.list=tau.list,gamma.list=gamma.list,n.fold=n.fold,pen.fac=pen.factors)
  tau <- res$tau
  gamma <- res$gamma
  b <- tlpreg0(y=Y,X=X.cbind,tau=tau,gamma=gamma,pen.fac=pen.factors)
  b1 <- b[1:(q1+q2)]
  b2 <- b[(q1+q2+1):p]
  nonzero.max <- min(sum(b2!=0),p-q1-q2-1)
  ord.b2 <- order(abs(b2),decreasing = TRUE)
  aic.memo <- numeric(nonzero.max+1)
  ind1 <- which(b1!=0)
  mod <- lm(Y~X.cbind[,ind1])
  aic.memo[1] <- AIC(mod)
  if (nonzero.max >0){
    for (i in 1:nonzero.max){
      ind2 <- which(abs(b2) >= abs(b2[ord.b2[i]]))
      ind <- c(ind1,ind2+q1+q2)
      mod <- lm(Y~X.cbind[,ind,drop=FALSE])
      aic.memo[i+1] <- AIC(mod)
    }
  }
  i <- which.min(aic.memo)-1
  if (i > 0){
    ind2 <- which(abs(b2) >= abs(b2[ord.b2[i]]))
    ind <- c(ind1,ind2+q1+q2)
  }else{
    ind <- ind1
  }
  b <- numeric(p)
  mod <- lm(Y~X.cbind[,ind])
  b[ind] <- as.numeric(mod$coefficients)[-1]
  return(list(b=b[1:q1],tau=tau,gamma=gamma))
}

coef.single.penalized <- function(Y,X,Z1,Z2,tau.list,gamma.list,n.fold){
  X.cbind <- cbind(X,Z1,Z2)
  p <- ncol(X.cbind)
  q1 <- ncol(X)
  q2 <- ncol(Z1)
  pen.factors <- rep(1,p)
  res <- cv.tlpreg0.aic(y=Y,X=X.cbind,tau.list=tau.list,gamma.list=gamma.list,n.fold=n.fold,pen.fac=pen.factors)
  tau <- res$tau
  gamma <- res$gamma
  b <- tlpreg0(y=Y,X=X.cbind,tau=tau,gamma=gamma,pen.fac=pen.factors)
  return(list(b=b[1:q1],tau=tau,gamma=gamma))
}

## This function finds relationships of parents, ancestors, interventions and 
## interventions of ancestors to be used in coefficient estimation
relations.finder <- function(U,W){
  p <- nrow(U)
  q <- nrow(W)
  
  U.depth <- topo.depth(U)
  
  U.parents <- t(U)
  U.ances <- matrix(0,p,p)
  U.int <- t(W)
  U.ancesint <- matrix(0,p,q)
  
  for (i in 0:max(U.depth)){
    idx.i <- which(U.depth==i)
    for (j in 1:length(idx.i)){
      idx <- idx.i[j]
      ances.mat <- U.ances[which(U.parents[idx,]==1),,drop=FALSE]
      U.ances[idx,] <- 1*((colSums(ances.mat)+U.parents[idx,,drop=FALSE])>0)
      parents.int <- U.int[which(U.parents[idx,]==1),,drop=FALSE]
      parents.ances <- U.ancesint[which(U.parents[idx,]==1),,drop=FALSE]
      U.ancesint[idx,] <- 1*(colSums(parents.int)+colSums(parents.ances) >0)
    }
  }
  return(list(U.parents,U.ances,U.int,U.ancesint))
}

## This function finds topological depth of primary variables,U should be a matrix
## with 0 representing zero elements and 1 representing non-zero elements
topo.depth <- function(U){
  p <- nrow(U)
  depth <- rep(p,p)
  dep.counter <- 0
  remained.idx <- (depth==p)
  while(sum(U)>0){
    if(min(colSums(U[remained.idx,remained.idx,drop=FALSE]))>0){
      stop(paste("Not DAG!"))
    }
    idx <- (colSums(U) == 0)
    dep <- rep(dep.counter,sum(idx))
    depth[idx] <- (depth[idx]+dep)/2 - abs(depth[idx]-dep)/2
    U[idx,] <- 0
    dep.counter <- dep.counter+1
    remained.idx <- (depth == p)
  }
  depth[depth==p] <- dep.counter
  return(depth)
}

## These two functions are to determine estimation order of elements in U, U should be a matrix
## with 0 representing zero elements and 1 representing non-zero elements
es.order <- function(U,depth){
  p <- nrow(U)
  order.mat <- outer(depth,depth,function(x,y){y*(y-1)/2 + y-x})
  order.mat <- order.mat*U
  order.unique <- unique(as.vector(order.mat))
  order.unique <- order.unique[order(order.unique)]
  order.unique.len <- length(order.unique)
  from.order.plus1 <- order.unique + 1
  to.order <- 0:(order.unique.len-1)
  convert.tab <- numeric(max(from.order.plus1))
  convert.tab[from.order.plus1] <- to.order
  order.mat <- apply(order.mat,1:2,function(x) {convert.tab[x+1]})
  return(order.mat)
}

es.order.sub <- function(order.mat){
  layer <- max(order.mat)
  p <- ncol(order.mat)
  idx.mat <- NULL
  for(i in 1:layer){
    idx <- which(order.mat==i)
    idx.col <- ceiling(idx/p)
    idx.row <- idx - (idx.col-1)*p
    submat <- cbind(rep(i,length(idx)),idx.row,idx.col)
    idx.mat <- rbind(idx.mat,submat)
  }
  return(idx.mat)
}
