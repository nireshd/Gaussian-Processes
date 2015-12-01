# Set of functions for performing Gaussian process regression.

############################################################################################################

# Calculate squared exponential covariance matrix.
# X1 and X2 - N1 x d and N2 x d data frames of the N1 and N2 input points in d dimensions that we wish to find 
# the covariance between.
# log.theta - vector of log parameters, i.e. (l1, ..., ld, sf, sn) where l1 to ld are the d length parameters, sf is 
# the signal variance, and sn is the noise variance.
sq.exp.cov <- function(X1,X2,theta){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  if(ncol(X1) != ncol(X2)){
    stop("X1 and X2 must contain input values of the same dimension.")
  } else if(!all(theta>=0)){
    stop("All parameters must be >= 0.")
  }
  N1 <- nrow(X1) # No. of points in X1.
  N2 <- nrow(X2) # No. of points in X2.
  d <- ncol(X1) # Dimension.
  sigma <-  matrix(rep(0, N1*N2),nrow=N1) # Matrix of zeroes,
  sf <- theta[d+1]
  sn <- theta[d+2]
  for(i in 1:d){
    l <- theta[i]
    xdiff <- (outer(X1[,i],X2[,i],function(x,y) x - y)/l)^2
    sigma <- sigma + xdiff
  }
  if(identical(X1,X2)){
    id <- diag(rep(1,N1))
    sigma.final <- sf*exp(-0.5*sigma) + sn*id
  } else {
    sigma.final <- sf*exp(-0.5*sigma)
  }
  return(sigma.final)
}

############################################################################################################

# Calculate inverse of the covariance function using svd.
cov.inverse.svd <- function(X, tol = sqrt(.Machine$double.eps)){
  # Generalized Inverse of a Matrix. Code from R function ginv.
  dnx <- dimnames(X)
  if(is.null(dnx)) dnx <- vector("list", 2)
  s <- svd(X)
  nz <- s$d > tol * s$d[1]
  K.inv <- structure(
    if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X,
    dimnames = dnx[2:1])
  # Logarithm of determinant.
  log.K.det <- sum(log(s$d))
  return(list(Inv = K.inv,lDet = log.K.det))
}

# Calculate inverse of the covariance function using Cholesky.
cov.inverse.chol <- function(X){
  R <- chol(X)
  Rt <- Conj(t(R))
  R.inv <- solve(R)
  Rt.inv <- solve(Rt)
  X.inv <- R.inv %*% Rt.inv
  log.X.det <- 2*sum(log(diag(R))) 
  return(list(Inv = X.inv, lDet = log.X.det))
}



############################################################################################################
# Functions for maximum likelihood estimation of parameters.

# Minus log marginal likelihood function of the Gaussian process with squared exponential covariance. 
# Used with optimisation to find parameter values.
# log.theta - vector of log parameters, i.e. (l1, ..., ld, sf, sn) where l1 to ld are the d length parameters, sf is 
# the signal variance, and sn is the noise variance.
# X - N x d data frame of N input points in d dimensions.
# Y - vector of N output values.
# inverse.method - "Cholesky" or "SVD".
mlml.sq.exp <- function(log.theta,X,Y, inverse.method="Cholesky"){
  theta <- exp(log.theta)
  K <- sq.exp.cov(X,X,theta)
  if(inverse.method=="Cholesky"){
    KK <- cov.inverse.chol(K)
  } else if(inverse.method=="SVD"){
    KK <- cov.inverse.svd(K)
  }
  mlml <- 0.5*(t(Y) %*% KK$Inv %*% Y) + 0.5*KK$lDet + 0.5*length(Y)*log(2*pi)
  return(as.numeric(mlml))
}

# Minus log marginal likelihood derivative function of the Gaussian process with squared exponential covariance. 
# Used with optimisation to find parameter values.
# log.theta - vector of log parameters, i.e. (l1, ..., ld, sf, sn) where l1 to ld are the d length parameters, sf is 
# the signal variance, and sn is the noise variance.
# X - N x d data frame of N input points in d dimensions.
# Y - vector of N output values.
# inverse.method - "Cholesky" or "SVD".
mlml.grad.sq.exp <- function(log.theta, X, Y, inverse.method="Cholesky"){
  X <- as.matrix(X)
  theta <- exp(log.theta)
  d <- ncol(X)
  K <- sq.exp.cov(X,X,theta)
  A <- sq.exp.cov(X,X,c(theta[1:d],1,0))
  if(inverse.method=="Cholesky"){
    KK <- cov.inverse.chol(K)
  } else if(inverse.method=="SVD"){
    KK <- cov.inverse.svd(K)
  }
  alpha <- (KK$Inv%*%Y)%*%t(KK$Inv%*%Y) - KK$Inv
  grad <- c()
  sf <- theta[d+1]
  sn <- theta[d+2]
  for(i in 1:d){
    l <- theta[i]
    xdiff <- (outer(X[,i],X[,i],function(x,y) x - y)^2)/(l*l*l)
    dl <- sf*A*xdiff
    dml <- -0.5*sum(diag(alpha%*%dl))*l
    grad <- c(grad,dml)
  }
  dml.sf <- -0.5*sum(diag(alpha%*%A))*sf
  dml.sn <- -0.5*sum(diag(alpha))*sn
  grad <- c(grad,dml.sf,dml.sn)
  return(grad)
}


############################################################################################################
# Functions for maximum posterior estimation of parameters where we have a Gamma prior on the last covariate length parameter.

# Minus log posterior likelihood function of the Gaussian process with squared exponential covariance and Gamma prior on the last length parameter. 
# Used with optimisation to find parameter values.
# log.theta - vector of log parameters, i.e. (l1, ..., ld, sf, sn) where l1 to ld are the d length parameters, sf is 
# the signal variance, and sn is the noise variance.
# X - N x d data frame of N input points in d dimensions.
# Y - vector of N output values.
# shape, scale - gamma prior parameter values.
# inverse.method - "Cholesky" or "SVD".
mlmp.sq.exp.gammalc <- function(log.theta,X,Y,shape,scale, inverse.method="Cholesky"){
  theta <- exp(log.theta)
  K <- sq.exp.cov(X,X,theta)
  if(inverse.method=="Cholesky"){
    KK <- cov.inverse.chol(K)
  } else if(inverse.method=="SVD"){
    KK <- cov.inverse.svd(K)
  }
  ls <- theta[(length(theta)-2)]
  mlml <- 0.5*(t(Y) %*% KK$Inv %*% Y) + 0.5*KK$lDet + 0.5*length(Y)*log(2*pi) - (shape-1)*log(ls) + (ls/scale) + shape*log(scale) + lgamma(shape)
  return(as.numeric(mlml))
}

# Minus log marginal likelihood derivative function of the Gaussian process with squared exponential covariance. 
# Used with optimisation to find parameter values.
# log.theta - vector of log parameters, i.e. (l1, ..., ld, sf, sn) where l1 to ld are the d length parameters, sf is 
# the signal variance, and sn is the noise variance.
# X - N x d data frame of N input points in d dimensions.
# Y - vector of N output values.
# shape, scale - gamma prior parameter values.
# inverse.method - "Cholesky" or "SVD".
mlmp.grad.sq.exp.gammalc <- function(log.theta, X, Y, shape, scale, inverse.method="Cholesky"){
  X <- as.matrix(X)
  theta <- exp(log.theta)
  d <- ncol(X)
  K <- sq.exp.cov(X,X,theta)
  A <- sq.exp.cov(X,X,c(theta[1:d],1,0))
  if(inverse.method=="Cholesky"){
    KK <- cov.inverse.chol(K)
  } else if(inverse.method=="SVD"){
    KK <- cov.inverse.svd(K)
  }
  alpha <- (KK$Inv%*%Y)%*%t(KK$Inv%*%Y) - KK$Inv
  grad <- c()
  sf <- theta[d+1]
  sn <- theta[d+2]
  for(i in 1:d){
    l <- theta[i]
    xdiff <- (outer(X[,i],X[,i],function(x,y) x - y)^2)/(l*l*l)
    dl <- sf*A*xdiff
    dml <- -0.5*sum(diag(alpha%*%dl))*l
    grad <- c(grad,dml)
  }
  grad[d] <- grad[d] - ((shape-1)/theta[d]) + (1/scale)
  dml.sf <- -0.5*sum(diag(alpha%*%A))*sf
  dml.sn <- -0.5*sum(diag(alpha))*sn
  grad <- c(grad,dml.sf,dml.sn)
  return(grad)
}


###############################################################################################################
# Function to perform prediction.
# Calculate predicted mean and covariance from Gaussian process regression given data and parameter values.
# Predictions for input Xstar, given data X with output Y, and parameters theta.
# inverse.method - "Cholesky" or "SVD".
gpr.prediction <- function(X.test, theta, X, Y, inverse.method="Cholesky", alpha=0.05){
  z <- qnorm(1 - alpha/2)
  Xstar <- as.matrix(X.test)
  d <- ncol(as.matrix(X))
  K <- sq.exp.cov(X,X,theta)
  if(inverse.method=="Cholesky"){
    KK <- cov.inverse.chol(K)
  } else if(inverse.method=="SVD"){
    KK <- cov.inverse.svd(K)
  }
  Ky <- KK$Inv %*% Y
  result <- apply(Xstar, 1, function(x){
    XX <- matrix(x,nrow=1)
    Kstar <- sq.exp.cov(X, XX, c(theta[1:(d+1)],0))
    Kstarstar <- sq.exp.cov(XX,XX,c(theta[1:(d+1)],0))
    mu <- t(Kstar) %*% Ky
    cv <- Kstarstar - (t(Kstar) %*% KK$Inv %*% Kstar)
    return(c(mu, mu - z*sqrt(cv), mu + z*sqrt(cv)))
  })
  prediction <- data.frame(t(result))
  colnames(prediction) <- c("Outcome", "Lower Confidence Interval","Upper Confidence Interval")
  return(prediction)
}
