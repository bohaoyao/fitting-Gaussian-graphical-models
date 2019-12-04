# since Y comes before W in the topological ordering, we have swapped the last two rows and columns in the correlation matrix
C <- matrix(data=c(1,0.516,0.453,0.322,0.332,0.516,1,0.438,0.405,0.417,0.453,0.438,1,0.596,0.538,0.322,0.405,0.596,1,0.541,0.332,0.417,0.538,0.541,1), nrow=5)
suffStat <- list(C=C, n=20700)
indepTest <- gaussCItest
summary(fci(suffStat, indepTest, alpha = 0.01, p=5))

#############################################################################################
# Simple DAG
#############################################################################################

LambdaA1 <- matrix(data = c(0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,1,0,0,1,1,1,1,0), nrow=5)
fitGaussianGraph(S=C, n=20700, LambdaA = LambdaA1, OmegaA = diag(5), learning_rate = 10^-3, convergence = 10^-10)

# Testing one fewer edge

summary(fci(suffStat, indepTest, alpha = 0.0001, p=5))
LambdaA2 <- matrix(data = c(0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,1,1,1,0), nrow=5)
Sigma1 <- fitGaussianGraph(S=C, n=20700, LambdaA = LambdaA1, OmegaA = diag(5), learning_rate = 10^-3, convergence = 10^-10)$Sigmahat
Sigma2 <- fitGaussianGraph(S=C, n=20700, LambdaA = LambdaA2, OmegaA = diag(5), learning_rate = 10^-3, convergence = 10^-10)$Sigmahat
s <- 2 * (loglikelihoodS(Sigma1, C, 20700) - loglikelihoodS(Sigma2, C, 20700))
pchisq(s, df=1, lower.tail = FALSE)

#############################################################################################
# ADMG - last two edges are bidirected
#############################################################################################

LambdaA3 <- matrix(data = c(0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,1,0,0,1,1,1,0,0), nrow=5)
OmegaA3 <- matrix(data = c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,1,1), nrow=5)
fitGaussianGraph(S=C, n=20700, LambdaA = LambdaA3, OmegaA = OmegaA3, learning_rate = 10^-3, convergence = 10^-10)

# Testing one fewer edge

LambdaA4 <- matrix(data = c(0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,1,1,0,0), nrow=5)
Sigma1 <- fitGaussianGraph(S=C, n=20700, LambdaA = LambdaA3, OmegaA = OmegaA3, learning_rate = 10^-3, convergence = 10^-10)$Sigmahat
Sigma2 <- fitGaussianGraph(S=C, n=20700, LambdaA = LambdaA4, OmegaA = OmegaA3, learning_rate = 10^-3, convergence = 10^-10)$Sigmahat
s <- 2 * (loglikelihoodS(Sigma1, C, 20700) - loglikelihoodS(Sigma2, C, 20700))
pchisq(s, df=1, lower.tail = FALSE)

#############################################################################################
# ADMG - father's edges become bidirected
#############################################################################################

LambdaA <- matrix(data = c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0), nrow=5)
OmegaA <- matrix(data = c(1,0,1,0,1,0,1,1,0,1,1,1,1,0,0,0,0,0,1,1,1,1,0,1,1), nrow=5)
fitGaussianGraph(S=C, n=20700, LambdaA = LambdaA, OmegaA = OmegaA, learning_rate = 10^-3, convergence = 10^-10)

# Testing one fewer edge

OmegaA2 <- matrix(data = c(1,0,1,0,0,0,1,1,0,1,1,1,1,0,0,0,0,0,1,1,0,1,0,1,1), nrow=5)
Sigma1 <- fitGaussianGraph(S=C, n=20700, LambdaA = LambdaA, OmegaA = OmegaA, learning_rate = 10^-3, convergence = 10^-10)$Sigmahat
Sigma2 <- fitGaussianGraph(S=C, n=20700, LambdaA = LambdaA, OmegaA = OmegaA2, learning_rate = 10^-3, convergence = 10^-10)$Sigmahat
s <- 2 * (loglikelihoodS(Sigma1, C, 20700) - loglikelihoodS(Sigma2, C, 20700))
pchisq(s, df=1, lower.tail = FALSE)
