#' @useDynLib grivet
#' @importFrom stats lm AIC

#' @title Matrix V Estimation Function
#' @description This function estimates matrix V.
#' @param X A n*q matrix for intervention variables.
#' @param Y A n*p matrix for primary variables.
#' @param taus Parameter used for TLP in each regression model where Y[,i] is the response with X the predictor matrix.
#' @param gammas Parameter used for TLP in each regression model where Y[,i] is the response with X the predictor matrix.
#' @return An estimate of matrix V.
#' @references Discovery and inference of a causal network with hidden confounding.
#' @examples
#' \dontrun{
#' to be added.
#' }
#' @export
intdag.pmle.diff.aic <- function(X,Y,taus,gammas){
  p <- ncol(Y)
  q <- ncol(X)
  V <- matrix(0,q,p)
  for (i in 1:p){
    y=Y[,i]
    b <- tlpreg0(y,X=X,tau=taus[i],gamma=gammas[i])
    ind <- b!=0
    b.ols <- matrix(0,q,1)
    if (sum(ind)==0){
      b.ols[,1] <- b
    }else{
      aic.memo <- numeric(sum(ind))
      ord <- order(abs(b),decreasing = TRUE)
      for (l in 1:sum(ind)){
        ind2 <- which(abs(b)>=abs(b[ord[l]]))
        mod <- lm(y~X[,ind2,drop=FALSE])
        aic.memo[l] <- AIC(mod)
      }
      l <- which.min(aic.memo)
      ind <- which(abs(b)>=abs(b[ord[l]]))
      mod <- lm(y~X[,ind,drop=FALSE])
      b.ols[ind,1] <- as.numeric(mod$coefficients)[-1]
    }
    V[,i] <- b.ols[,1]
  }
  return(V)
}

#' @title Matrix V Estimation Function with Cross-Validation
#' @description This function estimates matrix V with cross-validation.
#' @param X A n*q matrix for intervention variables.
#' @param Y A n*p matrix for primary variables.
#' @param tau.list A list of parameters to be selected for TLP in each regression model where Y[,i] is the response with X the predictor matrix.
#' @param gamma.list A list of parameters to be selected for TLP in each regression model where Y[,i] is the response with X the predictor matrix.
#' @param n.fold Number of folds to be used.
#' @return A list consists of an estimate of matrix V, taus and gammas selected
#' @references Discovery and inference of a causal network with hidden confounding.
#' @examples
#' \dontrun{
#' to be added.
#' }
#' @export
cv.intdag.pmle.diff.aic <- function(X,Y,tau.list,gamma.list,n.fold){
  p <- ncol(Y)
  taus <- numeric(p)
  gammas <- numeric(p)
  for (i in 1:p){
    res <- cv.tlpreg0.aic(Y[,i],X,tau.list,gamma.list,n.fold)
    taus[i] <- res$tau
    gammas[i] <- res$gamma
  }
  V = intdag.pmle.diff.aic(X,Y,taus,gammas)
  return(list(V=V,taus=taus,gammas=gammas))
}

#' @title Peeling Algorithm Implementation
#' @description This function conducts the peeling algorithm.
#' @param v A q*p matrix encoding relationships among Y and X.
#' @return A list consists of the estimates of ancestral relationships, interventional relationships and candidate-iv relationships.
#' @references Discovery and inference of a causal network with hidden confounding.
#' @examples
#' \dontrun{
#' to be added.
#' }
#' @export
topological_order <- function(v) {
  p <- ncol(v)
  q <- nrow(v)
  if (q < p) stop("No sufficient interventions: q < p.")
  
  an_mat <- matrix(0, p, p)
  in_mat <- matrix(0, q, p)
  iv_mat <- matrix(0, q, p)
  
  removed_x <- rep(FALSE, q)
  removed_y <- rep(FALSE, p)
  
  v_abs <- abs(v)
  v_nz <- v != 0
  
  # check if there is a primary variable without intervention.
  while (!all(removed_y)) {
    
    # leaf-instrument pairs
    iv_targets <- rowSums(as.matrix(v_nz[, !removed_y]))
    one <- min(iv_targets[iv_targets > 0 & !removed_x])
    leaf_iv <- which(!removed_x & iv_targets == one)
    if (length(leaf_iv) == 0) break
    leaf <- rep(NA, length(leaf_iv))
    leaf.iter <- 1
    for (l in leaf_iv) {
      ## j <- which(v_abs[l,] == max(v_abs[l, !removed_y]))[1] 
      j <- which((v_abs[l,] == max(v_abs[l, !removed_y]))*(!removed_y)==1)[1]
      iv_mat[l, j] <- 1
      in_mat[l, j] <- 1
      leaf[leaf.iter] <- j
      leaf.iter <- leaf.iter + 1
    }
    leaf <- unique(leaf)
    
    # leaf-noninstrument pairs
    for (j in leaf) {
      leaf_interventions <- which(!removed_x & v_abs[, j] != 0)
      leaf_noniv <- setdiff(leaf_interventions, leaf_iv)
      in_mat[leaf_noniv, j] <- 1
    }
    
    # ancestral relations
    for (j in leaf) {
      j_instrument <- which(iv_mat[, j] != 0)
      j_descendant <- vector("list", length = length(j_instrument))
      for (l in 1:length(j_instrument)) {
        l2 <- j_instrument[l]
        j_descendant[[l]] <- which(removed_y & v_abs[l2, ] != 0)
      }
      j_descendant <- Reduce(intersect, j_descendant)
      an_mat[j, j_descendant] <- 1
    }
    
    # peeling-off
    removed_y[leaf] <- TRUE
    removed_x[leaf_iv] <- TRUE
  }
  
  # reconstruction of topological order
  an_mat_interm <- (solve(diag(p) - an_mat) != 0) - diag(p)
  if(is.acyclic(an_mat_interm)){
    an_mat <- an_mat_interm
  }else{
    an_mat_interm <- an_mat
    for(l in 1:p){
      an_mat_interm <- an_mat_interm %*% an_mat
      an_mat <- an_mat + an_mat_interm
    }
  }
  in_mat <- 1*(in_mat%*%(diag(p)+an_mat)>0)
  
  list(an_mat = an_mat, in_mat = in_mat, iv_mat = iv_mat)
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