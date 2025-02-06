library(RANN)

# soFun.test -------------------------------------------------------------------------
#' Scalar on Function test for dependence.
#'
#' The test is based on the asymptotic distribution of the soFun statistic.
#' @param Y Vector (length n)
#' @param X Matrix (p by n)
#' @param eps Error bound for nearest neighbour search: a value of 0.0 implies exact nearest neighbour search.
#' @param na.rm Remove NAs if TRUE
#' @param permutation If TRUE, the p-value is calculated by permutation. If FALSE, the p-value is calculated based on the asymptotic distribution of the test statistic. The default is FALSE.
#' @param R Number of permutations. Ignored if permutation == FALSE.
#' @details The test statistic is the \code{\link{soFun}} statistic, with an appropriate scaling applied if permuatation == FALSE.
#' The test is based either on the known asymptotic normal distribution of this statistic or on a resampling procedure, according to the value of the argument \code{permutation}.
#' @return A list containing the test statistic and its p-value.
#' @import RANN
#' @export
#' @author Daniel Strenger
#' @references Hörmann, S. and Strenger, D. (2025). Azadkia–Chatterjee’s dependence coefficient for infinite dimensional data.
#' @examples
#' n <- 100
#' p <- 200
#' x <- brownian.motion(n,p)
#' y <- sin(2*pi* apply(x,2,mean))
#' if(soFun.test(y,x)$p.value < 0.05){
#'   print("The hypothesis of independence can be rejected at a significance level of 0.05.")
#' } else{
#'   print("The hypothesis of independence cannot be rejected at a significance level of 0.05.")
#' }
soFun.test <- function(Y, X, eps=2, na.rm=TRUE, permutation = FALSE, R = 1000){
  if(permutation == TRUE){
    statistic <- soFun(Y,X,eps,na.rm)
    statitstics.permuted <- numeric(R)
    for(i in 1:R){
      statitstics.permuted[i] <- soFun(sample(Y),X,eps,na.rm)
    }
    result <- list()
    result$statistic <- statistic
    result$p.value <- mean(statitstics.permuted > statistic)
    return(result)
  }

  if (na.rm == TRUE) {
    # NAs are removed here:
    ok <- complete.cases(Y,t(X))
    X <- as.matrix(X[,ok])
    Y <- Y[ok]
  }
  if(!is.vector(Y)) {
    Y <- as.vector(Y)
  }
  if(!is.matrix(X)) {
    X <- as.matrix(X)
  }

  n <- length(Y)

  Q <- soFun(Y,X)/6

  neighbor.indices <- nn2(t(X),k=2,eps=eps)$nn.idx[,2]

  #numbers of pairs of summands with 2,1, or 0 equal indices
  a2 <- sum(neighbor.indices[neighbor.indices]==(1:n))+n
  a1 <- sum(table(neighbor.indices)^2)+n-2*sum(neighbor.indices[neighbor.indices]==(1:n))
  a0 <- n^2-a1-a2

  #covariances of summands with 2,1, or 0 equal indices
  v1 <- (4*n^5-25*n^4+30*n^3+25*n^2-34*n)/(180*n*(n-1)*(n-2))
  v2 <- (n^3+-2*n^2-n+2)/(18*(n-1))
  v0 <- (-4*(n+1))/45

  #conditional variance of the test statistic given X before normalising
  v <- (a1*v1+a2*v2+a0*v0)/n^3

  result <- list()
  result$statistic <- sqrt(n)*Q/sqrt(v)
  result$p.value <- 1-pnorm(result$statistic)

  return(result)
}

# soFun -------------------------------------------------------------------------
#' Estimate the Scalar on Function Dependence Coefficient (soFun).
#'
#'The Scalar on Function Dependence Coefficient (soFun) is an adaption of the conditional dependence coefficient (CODEC) by Azadkia and Chatterjee to the setting of continuously distributed response Y and functional covariate X. It is a measure of the amount of dependence of Y on X, based on an i.i.d. sample of (Y, X). The coefficient is asymptotically guaranteed to be between 0 and 1.
#' @references Hörmann, S. and Strenger, D. (2024). Measuring dependence between a scalar response and a functional covariate.
#' @references Azadkia, M. and Chatterjee, S. (2019). A simple measure of conditional dependence. \href{https://arxiv.org/pdf/1910.12327.pdf}{https://arxiv.org/pdf/1910.12327.pdf}.
#' @details
#' The returned statistic can be positive or negative. Asymptotically, it is guaranteed to be between 0 and 1. A small value indicates low dependence of Y on X, and a high value indicates strong dependence. The statistic is the base for the \code{\link{soFun.test}}.
#' @export
#' @param X Matrix of functional covariates (p by n).
#' @param Y Vector of scalar responses (length n).
#' @param eps Error bound for nearest neighbour search: a value of 0.0 implies exact nearest neighbour search.
#' @param na.rm Remove NAs if TRUE
#' @return The measure of dependence of Y on X.
#' @examples
#' n <- 100
#' p <- 200
#' x <- brownian.motion(n,p)
#' y <- sin(2*pi* apply(x,2,mean))
#' soFun(y,x)

soFun <- function(Y, X, eps=2, na.rm=TRUE){
  if (na.rm == TRUE) {
    # NAs are removed here:
    ok <- complete.cases(Y,t(X))
    X <- as.matrix(X[,ok])
    Y <- Y[ok]
  }
  n <- length(Y)
  if(dim(X)[2]!=n){
    stop("Lengths of X and Y do not match!")
  }

  u <- rank(Y)/n

  neighbor.indices <- nn2(t(X),k=2,eps=eps)$nn.idx[,2]

  f <- cbind(u,u[neighbor.indices])
  f <- apply(f,1,min)
  return(6*mean(f)-2*(n+1)/n)
}
