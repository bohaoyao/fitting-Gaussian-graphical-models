###################
# Example 4.5
###################

Lambda <- matrix(c(0,0.5,0.25,0,0,0,0.5,0,0,0,0,1,0,0,0,0), nrow = 4, ncol = 4, byrow = TRUE)
Omega <- matrix(c(1,0,0,0,0,1,0,0.25,0,0,1,0,0,0.25,0,1), nrow = 4, ncol = 4, byrow = TRUE)
Sigma <- t(solve(diag(4)-Lambda)) %*% Omega %*% solve(diag(4)-Lambda)
U <- chol(Sigma)

n <- 1000 # number of samples
N <- 1000 # number of iterations
count <- 0

bootstrap1 <- function(U, n){
  Z <- matrix(rnorm(4*n),n,4)
  X <- Z %*% U
  S <- cov(X)
  VermaVector <- numeric(1000)
  for (i in 1:1000){
    v <- sample(1:n, n, replace = TRUE)
    X1 <- X[v,]
    S <- cov(X1) #sample covariance
    fVerma <- S[1,1]*S[1,3]*S[2,2]*S[3,4]-S[1,2]^2*S[1,3]*S[3,4]-
      S[1,1]*S[1,4]*S[2,2]*S[3,3]+S[1,2]^2*S[1,4]*S[3,3]-
      S[1,1]*S[1,3]*S[2,3]*S[2,4]+S[1,1]*S[1,4]*S[2,3]^2+
      S[1,2]*S[1,3]^2*S[2,4]-S[1,2]*S[1,3]*S[1,4]*S[2,3]
    VermaVector[i] <- fVerma
  }
  return(VermaVector)
}

set.seed(1)

VermaVector <- bootstrap1(U, n)
print(quantile(VermaVector, probs = c(0.025,0.975))) # Confidence Interval

# Repeat the bootstrap 1000 times

for (i in 1:N) {
  VermaVector <- bootstrap1(U, n)
  if (0<unname(quantile(VermaVector, .975)) & unname(quantile(VermaVector, .025)<0)) {
      count = count + 1}
}
print(count) #number of times 0 is in the confidence interval

###################
# Non-example
###################

Lambda <- matrix(c(0,0.5,0.25,0.2,0,0,0.5,0,0,0,0,1,0,0,0,0), nrow = 4, ncol = 4, byrow = TRUE)
Sigma <- t(solve(diag(4)-Lambda)) %*% Omega %*% solve(diag(4)-Lambda)
U <- chol(Sigma)
count <- 0

set.seed(1)

VermaVector <- bootstrap1(U, n)
print(quantile(VermaVector, probs = c(0.025,0.975))) # Confidence Interval

# Repeat the bootstrap 1000 times

for (i in 1:N) {
  VermaVector <- bootstrap1(U, n)
  if (0<unname(quantile(VermaVector, .975)) & unname(quantile(VermaVector, .025)<0)) {
    count = count + 1}
}
print(count) #number of times 0 is in the confidence interval
