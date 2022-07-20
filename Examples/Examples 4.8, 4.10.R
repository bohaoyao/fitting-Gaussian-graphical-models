##################
# Example 4.8
##################

library(BCD)

Lambda <- matrix(c(0,0.5,0.25,0,0,0,0.5,0,0,0,0,1,0,0,0,0), nrow = 4, ncol = 4, byrow = TRUE)
Omega <- matrix(c(1,0,0,0,0,1,0,0.25,0,0,1,0,0,0.25,0,1), nrow = 4, ncol = 4, byrow = TRUE)
Sigma <- t(solve(diag(4)-Lambda)) %*% Omega %*% solve(diag(4)-Lambda)
U <- chol(Sigma)

LambdaA <- matrix(c(0,1,1,0,0,0,1,0,0,0,0,1,0,0,0,0), nrow = 4, ncol = 4, byrow = TRUE)
OmegaA <- matrix(c(1,0,0,0,0,1,0,1,0,0,1,0,0,1,0,1), nrow = 4, ncol = 4, byrow = TRUE)

set.seed(1)

n <- 1000
N <- 1000
chisq <- c()
chisq2 <- c()

for (i in 1:N) {
  Z <- matrix(rnorm(4*n),n,4)
  X <- Z %*% U
  chisq[i] <- fitGaussianGraph(X=X, LambdaA=LambdaA, OmegaA=OmegaA, learning_rate = 10^-1, convergence_rate = 10^-10)$chisq
  out <- ricf(t(LambdaA), Omega, t(X))
  S <- t(X) %*% X/n
  chisq2[i] <- 2 * (loglikelihoodS(S, S, n) - loglikelihoodS(out$SigmaHat, S, n)) #comparism with RICF
  print(i)
}

ks.test(chisq, pchisq, df=1)
ks.test(chisq2, pchisq, df=1) #comparism with RICF
chisq1 <- ecdf(chisq)
Ex <- ecdf(qchisq(seq(0,0.9999,0.0001), df=1))
plot(chisq1, verticals=TRUE, do.points=FALSE)
plot(Ex, verticals=TRUE, do.points=FALSE, col='red', add=TRUE)


###########################################################
# re-run for small n (Example 4.10)
###########################################################

# n = 20

set.seed(1)

n <- 20
N <- 1000
chisq <- c()

# Warning: this code will take a while due to small learning rate. At a higher learning rate we may overshoot to a singular matrix.
for (i in 1:N) {
  Z <- matrix(rnorm(4*n),n,4)
  X <- Z %*% U
  chisq[i] <- fitGaussianGraph(X=X, LambdaA=LambdaA, OmegaA=OmegaA, learning_rate = 10^-3, convergence_rate = 10^-10)$chisq
  print(i)
}

ks.test(chisq, pchisq, df=1)
chisq1 <- ecdf(chisq)
Ex <- ecdf(qchisq(seq(0,0.9999,0.0001), df=1))
plot(chisq1, verticals=TRUE, do.points=FALSE)
plot(Ex, verticals=TRUE, do.points=FALSE, col='red', add=TRUE)

# re-run for RICF

set.seed(1)

n <- 20
N <- 1000
chisq <- c()

for (i in 1:N) {
  Z <- matrix(rnorm(4*n),n,4)
  X <- Z %*% U
  out <- ricf(t(LambdaA), Omega, t(X))
  S <- t(X) %*% X/n
  chisq[i] <- 2 * (loglikelihoodS(S, S, n) - loglikelihoodS(out$SigmaHat, S, n)) 
}

ks.test(chisq, pchisq, df=1)
chisq1 <- ecdf(chisq)
Ex <- ecdf(qchisq(seq(0,0.9999,0.0001), df=1))
plot(chisq1, verticals=TRUE, do.points=FALSE)
plot(Ex, verticals=TRUE, do.points=FALSE, col='red', add=TRUE)

# n = 10

set.seed(1)

n <- 10
N <- 1000
chisq <- c()

# Warning: this code will take a while due to small learning rate. At a higher learning rate we may overshoot to a singular matrix.
for (i in 1:N) {
  Z <- matrix(rnorm(4*n),n,4)
  X <- Z %*% U
  chisq[i] <- fitGaussianGraph(X=X, LambdaA=LambdaA, OmegaA=OmegaA, learning_rate = 10^-3, convergence_rate = 10^-10)$chisq
  print(i)
}

ks.test(chisq, pchisq, df=1)
chisq1 <- ecdf(chisq)
Ex <- ecdf(qchisq(seq(0,0.9999,0.0001), df=1))
plot(chisq1, verticals=TRUE, do.points=FALSE)
plot(Ex, verticals=TRUE, do.points=FALSE, col='red', add=TRUE)

# re-run for RICF

set.seed(1)

n <- 10
N <- 1000
chisq <- c()

for (i in 1:N) {
  Z <- matrix(rnorm(4*n),n,4)
  X <- Z %*% U
  out <- ricf(t(LambdaA), Omega, t(X))
  S <- t(X) %*% X/n
  chisq[i] <- 2 * (loglikelihoodS(S, S, n) - loglikelihoodS(out$SigmaHat, S, n)) 
}

ks.test(chisq, pchisq, df=1)
chisq1 <- ecdf(chisq)
Ex <- ecdf(qchisq(seq(0,0.9999,0.0001), df=1))
plot(chisq1, verticals=TRUE, do.points=FALSE)
plot(Ex, verticals=TRUE, do.points=FALSE, col='red', add=TRUE)

#comparison tests

llhd <- function(Sigma, S, n) {
  - n/2*(log(det(Sigma)) + sum(diag(solve.default(Sigma, S))))
}

set.seed(1)
n <- 10
for (i in 1:20) {
  Z <- matrix(rnorm(4*n),n,4)
  X <- Z %*% U
  S1 <- fitGaussianGraph(X=X, LambdaA=LambdaA, OmegaA=OmegaA, learning_rate = 10^-1, convergence_rate = 10^-12)$Sigmahat
  S2 <- ricf(t(LambdaA), OmegaA, t(X))$SigmaHat
  S <- t(X) %*% X/n
  print(llhd(S1, S, n)-llhd(S2, S, n))
}

set.seed(1)
n <- 1000
for (i in 1:20) {
  Z <- matrix(rnorm(4*n),n,4)
  X <- Z %*% U
  S1 <- fitGaussianGraph(X=X, LambdaA=LambdaA, OmegaA=OmegaA, learning_rate = 10^-1, convergence_rate = 10^-12)$Sigmahat
  S2 <- ricf(t(LambdaA), OmegaA, t(X))$SigmaHat
  S <- t(X) %*% X/n
  print(llhd(S1, S, n)-llhd(S2, S, n))
}
