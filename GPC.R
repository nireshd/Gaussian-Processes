# Set of functions for performing Gaussian process classification.
source("GPR.R")

######################################################################################
# Calculates the marginal log likelihood for expectation propagation given
# K - covariance matrix of the data
# sigma - approximate posterior covariance matrix
# mu - approximate posterior mean
# L - cholesky decomposition of B = I + S^(1/2) K S^(1/2)
# tau.tilde - inverse site variances
# nu.tilde - diag(tau.tilde) x site means
mll.ep <- function(K, sigma, mu, L, tau.tilde, nu.tilde){
  tau.cavity <- (1/diag(sigma)) - tau.tilde # Compute approximate cavity parameters.
  nu.cavity <- (mu/diag(sigma)) - nu.tilde
  S <- diag(tau.tilde) 
  t <- diag(tau.cavity)
  V <- solve(t(L))%*%(sqrt(S)%*%K)
  mu.cavity <- nu.cavity/tau.cavity
  term3 <- sum(pnorm((y*mu.cavity)/(sqrt(1 + (1/tau.cavity))),log.p=TRUE)) # Third term
  term1.4 <- 0.5*sum(log(1 + (tau.tilde/tau.cavity))) - sum(log(diag(L))) # First and fourth term
  term2.5a <- 0.5*t(nu.tilde)%*%((K - (t(V)%*%V) - solve(t + S))%*%nu.tilde) # Second and part of fifth term
  term5b <- 0.5*t(mu.cavity)%*%(t%*%(solve(S + t)%*%((S%*%mu.cavity)-(2*nu.tilde)))) # Rest of firth term
  mll <- term1.4 + term2.5a + term3 + term5b # Combine to make log marginal likelihood
  return(mll)
}

#######################################################################################
# Apply expectation propagation to data.
# Returns natural site parameters and marginal log likelihood.
# K - covariance matrix of data.
# y - vector of data classifications (either +1 or -1).
# tol - tolerance for convergence.
# max.sweeps - maximum number of sweeps.
expectation.propagation <- function(K, y, tol=sqrt(.Machine$double.eps), max.sweeps=10){
  n <- length(y)
  I <- diag(rep(1,n))
  nu.tilde <- rep(0,n) # Initialise parameters.
  tau.tilde <- rep(0,n)
  sigma <- K
  mu <- rep(0,n)
  lZ <- mll.ep(K, sigma, mu, diag(rep(1,n)), tau.tilde, nu.tilde) # Calculate initial marginal log likelihood.
  lZ.old <- Inf
  sweep <- 0
  conv <- "Not converged"
  while((abs(lZ - lZ.old) > tol) & (sweep < max.sweeps)){ # Repeat until log likelihood converges or reach max number of sweeps.
    lZ.old <- lZ
    sweep <- sweep + 1
    for(i in 1:n){ # Loop over data.
      tau.cavity <- (1/sigma[i,i]) - tau.tilde[i] # Compute approximate cavity parameters.
      nu.cavity <- (mu[i]/sigma[i,i]) - nu.tilde[i]
      z <- (y[i]*(nu.cavity/tau.cavity))/sqrt(1 + (1/tau.cavity)) # Compute marginal moments.
      mu.hat <- (nu.cavity/tau.cavity) + 
        ((y[i]*(1/tau.cavity)*dnorm(z))/(pnorm(z)*sqrt(1 + (1/tau.cavity))))
      sigma.hat <- (1/tau.cavity) - 
        (((1/tau.cavity^2)*dnorm(z))/((1 + (1/tau.cavity))*pnorm(z)))*(z + (dnorm(z)/pnorm(z)))
      delta.tau <- (1/sigma.hat) - tau.cavity - tau.tilde[i] # Update site parameters
      tau.tilde[i] <- tau.tilde[i] + delta.tau
      nu.tilde[i] <- (1/sigma.hat)*mu.hat - nu.cavity
      sigma <- sigma - (1/((1/delta.tau) + sigma[i,i]))*(sigma[,i]%*%t(sigma[,i])) # Update mu and sigma
      mu <- sigma %*% nu.tilde
    }
    S <- diag(tau.tilde)
    L <- chol(I + sqrt(S)%*%(K%*%sqrt(S))) # Re-compute approximate posterior parameters
    V <- solve(t(L))%*%(sqrt(S)%*%K)
    sigma <- K - t(V)%*%V
    mu <- sigma%*%nu.tilde
    lZ <- mll.ep(K, sigma, mu, L, tau.tilde, nu.tilde) # Re-compute marginal log likelihood.
    if(abs(lZ - lZ.old) <= tol){
      conv <- "Converged"
    }
  }
  print(conv)
  return(list(logZ=lZ, tau.site=tau.tilde, nu.site=nu.tilde, converged=conv))
}

#######################################################################################
# Make EP predictions from data with squared exponential covariance.
# Return vector of predictive class probabilites for class +1.
# X.test - test points to evaluate class probabilities.
# X - data inputs.
# Y - vector of data classifications (+1 or -1).
# theta - vector of d+1 hyperparameters for covariance function (d is input dimension).
prediction.ep.sq.exp <- function(X.test, X, Y, theta, eps=sqrt(.Machine$double.eps), max=10){
  n <- length(Y)
  K <- sq.exp.cov(X,X,c(theta,0))
  EP <- expectation.propagation(K, Y, tol=eps, max.sweeps=max)
  nu.tilde <- EP$nu.site
  tau.tilde <- EP$tau.site
  I <- diag(rep(1,n))
  S <- diag(tau.tilde)
  L <- chol(I + (sqrt(S)%*%(K%*%sqrt(S)))) # B = I + S^(1/2) K S^(1/2)
  z <- sqrt(S)%*%(solve(t(L))%*%(solve(L)%*%(sqrt(S)%*%(K%*%nu.tilde))))
  X.star <- as.matrix(X.test)
  pi.star <- apply(X.star, 1, function(x){
    XX <- matrix(x,nrow=1)
    Kstar <- sq.exp.cov(X, XX, c(theta,0))
    Kstarstar <- sq.exp.cov(XX,XX,c(theta,0))
    f.star <- t(Kstar)%*%(nu.tilde-z) # Calculate posterior mean.
    v <- solve(L)%*%(sqrt(S)%*%Kstar)
    v.star <- Kstarstar - t(v)%*%v # Calculate posterior variance.
    pi <- pnorm(f.star/sqrt(1 + v.star)) # Calculate predictive class probability for class +1.
    return(pi)
  })
  return(pi.star)
}

#######################################################################################
# Minus log marginal likelihood for squared exponential.
# log.theta - log of the covariance hyperparameters.
# X - data inputs.
# Y - vector of data classifications (+1 or -1).
# eps - tolerance for convergence of expectation propagation algorithm.
# max - maximum number of sweeps for expectation propagation algorithm.
mlml.ep.sq.exp <- function(log.theta, X, Y, eps=sqrt(.Machine$double.eps), max=10){
  theta <- exp(log.theta)
  K <- sq.exp.cov(X,X,c(theta,0)) # Calculate covariance.
  EP <- expectation.propagation(K, Y, tol=eps, max.sweeps=max) # Run EP algorithm.
  return(-EP$logZ)
}

# Minus log marginal likelihood gradient for squared exponential.
# log.theta - log of the covariance hyperparameters.
# X - data inputs.
# Y - vector of data classifications (+1 or -1).
# eps - tolerance for convergence of expectation propagation algorithm.
# max - maximum number of sweeps for expectation propagation algorithm.
mlml.grad.ep.sq.exp <- function(log.theta, X, Y, eps=sqrt(.Machine$double.eps), max=10){
  X <- as.matrix(X)
  n <- length(Y)
  d <- ncol(X)
  theta <- exp(log.theta)
  K <- sq.exp.cov(X,X,c(theta,0)) # Calculate covariance.
  EP <- expectation.propagation(K, Y, tol=eps, max.sweeps=max) # Run EP algorithm.
  S <- diag(EP$tau.site)
  I <- diag(rep(1,n))
  L <- chol(I + sqrt(S)%*%(K%*%sqrt(S))) # B = I + S^(1/2) K S^(1/2)
  b <- EP$nu.site - sqrt(S)%*%(solve(L)%*%(solve(t(L))%*%(sqrt(S)%*%(K%*%EP$nu.site))))
  R <- b%*%t(b) - sqrt(S)%*%(solve(t(L))%*%(solve(L)%*%sqrt(S))) # R = bb^T - S^(1/2) B^(-1) S^(1/2)
  sf <- theta[d+1]
  grad <- c()
  for(j in 1:d){ 
    l <- theta[j]
    xdiff <- (outer(X[,j],X[,j],function(x,y) x - y)^2)/(l*l*l)
    dl <- K*xdiff # Derivative matrix.
    dml <- -0.5*sum(diag(R%*%dl))*l # Derivative.
    grad <- c(grad,dml)
  }
  dml.sf <- -0.5*sum(diag(R%*%(K/sf)))*sf # Derivative.
  grad <- c(grad,dml.sf)
  return(grad)
}