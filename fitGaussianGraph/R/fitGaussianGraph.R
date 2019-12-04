#outputs value of a_{i,j}
f <- function(i, j, LambdaE, S){
  ans <- S[i,j]
  if (i > 1) {
    for (k in 1:(i-1)) {
      ans <- ans - LambdaE[k,i]*S[k,j]
    }
  }
  return(ans)
}

#outputs parents of i
parents <- function(i, LambdaA){
  ans <- c()
  if (i > 1) {
    for (j in 1:(i-1)) {
      if (LambdaA[j,i] != 0) {
        ans <- c(ans, j)
      }
    }
  }
  return(ans)
}

#outputs the starting point \Tilde{\Lambda}
SolveLambdaE <- function(LambdaA, S){
  p <- dim(S)[1] # dimension
  LambdaE <- matrix(0, p, p)
  for (i in 2:p) {
    pa <- parents(i, LambdaA)
    p1 <- length(pa)
    if (p1 > 0) {
      A <- matrix(0, p1, p1)
      b <- matrix(0, p1, 1)
      count1 <- 1
      for (j in pa) {
        count2 <- 1
        for (k in pa) {
          A[count1, count2] <- f(j, k, LambdaE, S)
          count2 <- count2 + 1
        }
        b[count1, 1] <- f(j, i, LambdaE, S)
        count1 <- count1 + 1
      }
      B <- solve(A,b)
      LambdaE[pa,i] <- B
    }
  }
  return(LambdaE)
}

#outputs the value of \Tilde{\omega_{ij}}
SolveOmega <- function(i, j, LambdaE, S){
  ans <- f(i, j, LambdaE, S)
  if (j > 1) {
    for (k in 1:(j-1)) {
      ans <- ans - LambdaE[k,j]*f(i, k, LambdaE, S)
    }
  }
  return(ans)
}

#outputs the starting points \Tilde{\Omega}
SolveOmegaE <- function(OmegaA, LambdaE, S){
  p <- dim(S)[1] # dimension
  OmegaE <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in i:p) {
      if (OmegaA[i,j] != 0) {
        OmegaE[i,j] <- SolveOmega(i, j, LambdaE, S)
        if (i != j) {
          OmegaE[j,i] <- OmegaE[i,j]
        }
      }
    }
  }
  return(OmegaE)
}

#outputs the derivative of the loglikelihood against Lambda
dLambda <- function(L, O, LambdaA, S) {
  p <- dim(O)[1] # dimension
  A <- 2*S%*%(diag(p)-L)%*%solve(O)
  for (i in 1:p) {
    for (j in 1:p) {
      if (LambdaA[i,j] == 0){
        A[i,j] <- 0
      }
    }
  }
  return(A)
}

#outputs the derivative of the loglikelihood against Omega
dOmega <- function(L, O, OmegaA, S) {
  p <- dim(O)[1] # dimension
  B <- solve(O)%*%t((diag(p)-L))%*%S%*%(diag(p)-L)%*%solve(O)-solve(O)
  for (i in 1:p) {
    for (j in i:p) {
      if (OmegaA[i,j] == 0){
        B[i,j] <- 0
        B[j,i] <- 0
      }
    }
  }
  return(B)
}

#outputs the loglikelihood given Lambda and Omega (divided by n)
loglikelihood <- function(L, O, S) {
  p <- dim(O)[1] # dimension
  return(log(det((diag(p)-L)%*%solve(O)%*%t(diag(p)-L)))-sum(diag(S%*%(diag(p)-L)%*%solve(O)%*%t(diag(p)-L))))
}

#outputs the loglikelihood given Sigma
loglikelihoodS <- function(Sigma, S, n) {
  return((-n*(log(det(Sigma)))-n*sum(diag(solve(Sigma)%*%S)))/2)
}

#ouputs the number of missing edges in the graph (used for finding the degree of freedom)
count_missing_edges <- function(LambdaA, OmegaA) {
  p <- dim(LambdaA)[1]
  edges <- 0
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      if (LambdaA[i,j] != 0) {
        edges <- edges + 1
      }
      else if (OmegaA[i,j] != 0) {
        edges <- edges + 1
      }
    }
  }
  return(p*(p-1)/2 - edges)
}

#uses hill climbing algorithm to output the supremum of Lambda, Omega and Sigma
#also computes the p-value
HillClimb <- function(LambdaE, OmegaE, learning_rate, convergence, S, n) {
  diff <- 1
  Lambda0 <- LambdaE
  Omega0 <- OmegaE
  p <- dim(S)[1] # dimension

  min_eig <- min(eigen(Omega0)$values)

  if (min_eig <= 0) {
    Omega0 <- Omega0 - (min_eig - 1)*diag(p)
  }


#  if (min(eigen(Omega0)$values) <= 0) {
#    SigmaE <- solve(t(diag(p)-LambdaE)) %*% OmegaE %*% solve(diag(p)-LambdaE)
#    message('Use a larger sample size to obtain p-value.')
#    return(list(Lambdahat = LambdaE, Omegahat = OmegaE, Sigmahat = SigmaE))
#  }

  likelihoodvector <- c(loglikelihood(Lambda0, Omega0, S))

  while (diff > convergence) {
    LambdaE <- Lambda0
    OmegaE <- Omega0
    A <- dLambda(LambdaE, OmegaE, LambdaE, S)
    B <- dOmega(LambdaE, OmegaE, OmegaE, S)
    Lambda0 <- LambdaE + learning_rate*A
    Omega0 <- OmegaE + learning_rate*B
    diff <- loglikelihood(Lambda0, Omega0, S) - loglikelihood(LambdaE, OmegaE, S)
    likelihoodvector <- c(likelihoodvector, loglikelihood(Lambda0, Omega0, S))
  }
  SigmaE <- solve(t(diag(p)-LambdaE)) %*% OmegaE %*% solve(diag(p)-LambdaE)
  s <- 2 * (loglikelihoodS(S, S, n) - loglikelihoodS(SigmaE, S, n)) #likelihood ratio statistic
  e <- count_missing_edges(LambdaE, OmegaE)
  plot(likelihoodvector)
  return(list(Lambdahat = LambdaE, Omegahat = OmegaE, Sigmahat = SigmaE, chisq = s, pvalue = pchisq(s, df=e, lower.tail = FALSE)))
}

#' Fitting bridgeless Gaussian graphical model
#'
#' Fits a bridgeless Gaussian graphical model as well as outputting the p-value of whether the data set fits the model.
#'
#' The p-value should be used as a guide for rejecting models only as we might no longer be working in a concave space at low p-values.
#' @param X centered data set (can provide both S and n instead).
#' @param S covariance matrix of X (optional if X is provided).
#' @param n sample size (optional if X is provided).
#' @param LambdaA adjacency matrix for the directed edges in G.
#' @param OmegaA adjacency matrix for the bidirected edges in G.
#' @param learning_rate learning rate for the hill-climbing algorithm.
#' @param convergence_rate the difference of the objective function between iterations whence we determine the hill-climbing algorithm has converged.
#'
#' @return A list with components
#' \describe{
#'   \item{Lambdahat}{a square matrix of the fitted regression coefficients.}
#'   \item{Omegahat}{the fitted covariance matrix of the error terms.}
#'   \item{Sigmahat}{the fitted covariance matrix of all observed variables.}
#'   \item{chisq}{the test statistic which converges to a chi-squared distribution.}
#'   \item{pvalue}{the p-value of whether the data set X fits the graphical model.}
#' }
#'
#' @export

fitGaussianGraph <- function(X, S, n, LambdaA, OmegaA, learning_rate, convergence_rate = 0) {
  if (missing(S)) {
    S <- cov(X)
  }
  if (missing(n)) {
    n <- length(X[,1])
  }
  LambdaE <- SolveLambdaE(LambdaA, S)
  OmegaE <- SolveOmegaE(OmegaA, LambdaE, S)
  return(HillClimb(LambdaE, OmegaE, learning_rate, convergence_rate, S, n))
}
