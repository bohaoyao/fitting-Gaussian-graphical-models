eps <- 10^-1
a <- (2-2*eps)/(1+sqrt(5))
Lambda <- matrix(c(0,0,0.25,0,0,0,0,0.25,0,0,0,0,0,0,0,0), nrow = 4, ncol = 4, byrow = TRUE)
Omega <- matrix(c(1,a,0,a,a,1,a,0,0,a,1,0,a,0,0,1), nrow = 4, ncol = 4, byrow = TRUE)
Sigma <- t(solve(diag(4)-Lambda)) %*% Omega %*% solve(diag(4)-Lambda)
U <- chol(Sigma)

LambdaA <- matrix(c(0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0), nrow = 4, ncol = 4, byrow = TRUE)
OmegaA <- matrix(c(1,1,0,1,1,1,1,0,0,1,1,0,1,0,0,1), nrow = 4, ncol = 4, byrow = TRUE)

set.seed(1)

n <- 1000
N <- 1000
chisq <- c()

for (i in 1:N) {
  Z <- matrix(rnorm(4*n),n,4)
  X <- Z %*% U
  chisq[i] <- fitGaussianGraph(X=X, LambdaA=LambdaA, OmegaA=OmegaA, learning_rate=10^-3, convergence_rate=10^-10)$chisq
  print(i)
}

ks.test(chisq, pchisq, df=1)
chisq1 <- ecdf(chisq)
Ex <- ecdf(qchisq(seq(0,0.9999,0.0001), df=1))
plot(chisq1, verticals=TRUE, do.points=FALSE)
plot(Ex, verticals=TRUE, do.points=FALSE, col='red', add=TRUE)

###########################################################################################

eps <- 10^-2
a <- (2-2*eps)/(1+sqrt(5))
Lambda <- matrix(c(0,0,0.25,0,0,0,0,0.25,0,0,0,0,0,0,0,0), nrow = 4, ncol = 4, byrow = TRUE)
Omega <- matrix(c(1,a,0,a,a,1,a,0,0,a,1,0,a,0,0,1), nrow = 4, ncol = 4, byrow = TRUE)
Sigma <- t(solve(diag(4)-Lambda)) %*% Omega %*% solve(diag(4)-Lambda)
U <- chol(Sigma)

LambdaA <- matrix(c(0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0), nrow = 4, ncol = 4, byrow = TRUE)
OmegaA <- matrix(c(1,1,0,1,1,1,1,0,0,1,1,0,1,0,0,1), nrow = 4, ncol = 4, byrow = TRUE)

set.seed(1)

n <- 1000000
N <- 1000
chisq <- c()

for (i in 1:N) {
  Z <- matrix(rnorm(4*n),n,4)
  X <- Z %*% U
  chisq[i] <- fitGaussianGraph(X=X, LambdaA=LambdaA, OmegaA=OmegaA, learning_rate=10^-5, convergence_rate=10^-10)$chisq
  print(i)
}

ks.test(chisq, pchisq, df=1)
chisq1 <- ecdf(chisq)
Ex <- ecdf(qchisq(seq(0,0.9999,0.0001), df=1))
plot(chisq1, verticals=TRUE, do.points=FALSE)
plot(Ex, verticals=TRUE, do.points=FALSE, col='red', add=TRUE)
