Sigma <- matrix(c(1,.8,.4375,.76,.8,1,.8,.5,.4375,.8,1,.1,.76,.5,.1,1), nrow = 4, ncol = 4, byrow = TRUE)
U <- chol(Sigma)

LambdaA <- matrix(c(0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0), nrow = 4, ncol = 4, byrow = TRUE)
OmegaA <- matrix(c(1,0,1,1,0,1,0,1,1,0,1,0,1,1,0,1), nrow = 4, ncol = 4, byrow = TRUE)

set.seed(1)

n <- 1000
N <- 1000
chisq <-  rep(0, n)

for (i in 1:N) {
  Z <- matrix(rnorm(4*n),n,4)
  X <- Z %*% U
  chisq[i] <- tryCatch(fitGaussianGraph(X=X, LambdaA=LambdaA, OmegaA=OmegaA, learning_rate=10^-3, convergence_rate=10^-10)$chisq,
          error=function(e)
          tryCatch(fitGaussianGraph(X=X, LambdaA=LambdaA, OmegaA=OmegaA, learning_rate=10^-20, convergence_rate=10^-10)$chisq,
            error = function(e) NA))
  print(i)
}

sum(is.na(chisq))

set.seed(1)

n <- 1000000
N <- 1000
chisq <- c()

for (i in 1:N) {
  Z <- matrix(rnorm(4*n),n,4)
  X <- Z %*% U
  chisq[i] <- tryCatch(fitGaussianGraph(X=X, LambdaA=LambdaA, OmegaA=OmegaA, learning_rate=10^-3, convergence_rate=10^-10)$chisq,
                       error=function(e)
                         tryCatch(fitGaussianGraph(X=X, LambdaA=LambdaA, OmegaA=OmegaA, learning_rate=10^-20, convergence_rate=10^-10)$chisq,
                                  error = function(e) NA))
  print(i)
}

sum(is.na(chisq))
