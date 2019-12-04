###################
# Example 4.6
###################

Lambda <- matrix(c(0,0.5,0.25,0,0,0,0.5,0,0,0,0,1,0,0,0,0), nrow = 4, ncol = 4, byrow = TRUE)
LambdaA <- matrix(c(0,1,1,1,0,0,1,0,0,0,0,1,0,0,0,0), nrow = 4, ncol = 4, byrow = TRUE) #assume there is an edge between 1 and 4.
Omega <- matrix(c(1,0,0,0,0,1,0,0.25,0,0,1,0,0,0.25,0,1), nrow = 4, ncol = 4, byrow = TRUE)
Sigma <- t(solve(diag(4)-Lambda)) %*% Omega %*% solve(diag(4)-Lambda)
U <- chol(Sigma)

n <- 1000 # number of samples
N <- 1000 # number of iterations
count <- 0

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

#outputs the estimation for Lambda
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

bootstrap2 <- function(U, n){
  Z <- matrix(rnorm(4*n),n,4)
  X <- Z %*% U
  S <- cov(X)
  VermaVector <- numeric(1000)
  for (i in 1:1000){
    v <- sample(1:n, n, replace = TRUE)
    X1 <- X[v,]
    S <- cov(X1) #sample covariance
    LambdaE <- SolveLambdaE(LambdaA, S)
    VermaVector[i] <- LambdaE[1,4]
  }
  return(VermaVector)
}

set.seed(1)

VermaVector <- bootstrap2(U, n)
print(quantile(VermaVector, probs = c(0.025,0.975))) # Confidence Interval

# Repeat the bootstrap 1000 times

for (i in 1:N) {
  VermaVector <- bootstrap2(U, n)
  if (0<unname(quantile(VermaVector, .975)) & unname(quantile(VermaVector, .025)<0)) {
    count = count + 1}
}
print(count) #number of times 0 is in the confidence interval
