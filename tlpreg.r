dyn.load("ctlpreg.so")
library(cvTools)

lasso0 <- function(y, X, b.init=rep(0,ncol(X)), lambda, pen.fac=rep(1,ncol(X)), tol=1e-5, cd.maxit=1e+4) {

    n <- nrow(X)
    p <- ncol(X)
    .C('classo', y = as.double(y),
                 X = as.double(X),
                 b0 = as.double(numeric(1)),
                 b = as.double(b.init),
                 r = as.double(y),
                 n = as.integer(n),
                 p = as.integer(p),
                 lambda = as.double(lambda),
                 pen_fac = as.integer(pen.fac),
                 tol = as.double(tol),
                 cd_maxit = as.integer(cd.maxit))$b
}

tlpreg0 <- function(y, X, b.init=rep(0,ncol(X)), tau=0.01, gamma=0.5, pen.fac=rep(1,ncol(X)), tol=1e-5, dc.maxit=20, cd.maxit=1e+4) {

    n <- nrow(X)
    p <- ncol(X)
    .C('ctlpreg0', y = as.double(y),
                   X = as.double(X),
                   b0 = as.double(numeric(1)),
                   b = as.double(b.init),
                   r = as.double(y),
                   n = as.integer(n),
                   p = as.integer(p),
                   tau = as.double(tau),
                   gamma = as.double(gamma),
                   pen_fac = as.integer(pen.fac),
                   tol = as.double(tol),
                   dc_maxit = as.integer(dc.maxit),
                   cd_maxit = as.integer(cd.maxit))$b
}

cv.tlpreg0 <- function(y,X,tau.list,gamma.list,n.fold,b.init=rep(0,ncol(X)),pen.fac=rep(1,ncol(X)),tol=1e-5,dc.maxit=20,cd.maxit=1e+4){
    tau.len <- length(tau.list)
    gamma.len <- length(gamma.list)
    cvm <- matrix(0,tau.len,gamma.len)
    n <- nrow(X)
    q <- ncol(X)
    for (i in 1:tau.len){
        tau <- tau.list[i]
        for (j in 1:gamma.len){
            gamma <- gamma.list[j]
            cvm.sub <- numeric(n.fold)
            folds <- cvFolds(n,n.fold)
            for (k in 1:n.fold){
                X.train <- X[folds$subsets[folds$which!=k],]
                y.train <- y[folds$subsets[folds$which!=k]]
                X.test <- X[folds$subsets[folds$which==k],]
                y.test <- y[folds$subsets[folds$which==k]]
                b.ols <- matrix(0,q,1)
                b <- tlpreg0(y=y.train,X=X.train,tau=tau,gamma=gamma,b.init=b.init,pen.fac=pen.fac,tol=tol,dc.maxit=dc.maxit,cd.maxit=cd.maxit)
                ind <- b!=0
                if (sum(ind)==0){
                    b.ols[,1] <- b
                }else{
                    mod <- lm(y.train~X.train[,ind,drop=FALSE])
                    b.ols[ind,1] <- as.numeric(mod$coefficients)[-1]
                }
                cvm.sub[k] <- mean((y.test-X.test%*%b.ols)^2)
            }
            cvm[i,j] <- mean(cvm.sub)
        }
    }
    locations <- as.numeric(which(cvm==min(cvm,na.rm = TRUE),arr.ind = TRUE)[1,])
    tau <- tau.list[locations[1]]
    gamma <- gamma.list[locations[2]]
    return(list(tau=tau,gamma=gamma))
}

cv.tlpreg0.aic <- function(y,X,tau.list,gamma.list,n.fold,b.init=rep(0,ncol(X)),pen.fac=rep(1,ncol(X)),tol=1e-5,dc.maxit=20,cd.maxit=1e+4){
    tau.len <- length(tau.list)
    gamma.len <- length(gamma.list)
    cvm <- matrix(0,tau.len,gamma.len)
    n <- nrow(X)
    q <- ncol(X)
    for (i in 1:tau.len){
        tau <- tau.list[i]
        for (j in 1:gamma.len){
            gamma <- gamma.list[j]
            cvm.sub <- numeric(n.fold)
            folds <- cvFolds(n,n.fold)
            for (k in 1:n.fold){
                X.train <- X[folds$subsets[folds$which!=k],,drop=FALSE]
                y.train <- y[folds$subsets[folds$which!=k]]
                X.test <- X[folds$subsets[folds$which==k],,drop=FALSE]
                y.test <- y[folds$subsets[folds$which==k]]
                b.ols <- matrix(0,q,1)
                b <- tlpreg0(y=y.train,X=X.train,tau=tau,gamma=gamma,b.init=b.init,pen.fac=pen.fac,tol=tol,dc.maxit=dc.maxit,cd.maxit=cd.maxit)
                ind <- b!=0
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
                    mod <- lm(y.train~X.train[,ind,drop=FALSE])
                    b.ols[ind,1] <- as.numeric(mod$coefficients)[-1]
                }
                cvm.sub[k] <- mean((y.test-X.test%*%b.ols)^2)
            }
            cvm[i,j] <- mean(cvm.sub)
        }
    }
    try(locations <- as.numeric(which(cvm==min(cvm,na.rm = TRUE),arr.ind = TRUE)[1,]),{warning("You may change the penalty parameters used");
        locations <- c(ceiling(length(tau.list)/2),ceiling(length(gamma.list)/2))})
    tau <- tau.list[locations[1]]
    gamma <- gamma.list[locations[2]]
    return(list(tau=tau,gamma=gamma))
}
