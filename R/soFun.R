library(RANN)
library(stats)

# soFun -------------------------------------------------------------------------
#' Estimate the Scalar on Function Dependence Coefficient (soFun).
#'
#'The Scalar on Function Dependence Coefficient (soFun) is an adaption of the conditional dependence coefficient (CODEC) by Azadkia and Chatterjee to the setting of continuously distributed response Y and functional covariate X. It is a measure of the amount of dependence of Y on X, based on an i.i.d. sample of (Y, X). The coefficient is asymptotically guaranteed to be between 0 and 1.
#' @references Hörmann, S. and Strenger, D. (2025). Azadkia–Chatterjee’s dependence coefficient for infinite dimensional data. \href{https://arxiv.org/abs/2405.07732}{https://arxiv.org/abs/2405.07732}
#' @references Azadkia, M. and Chatterjee, S. (2019). A simple measure of conditional dependence. \href{https://arxiv.org/pdf/1910.12327.pdf}{https://arxiv.org/pdf/1910.12327.pdf}.
#' @details
#' The returned statistic can be positive or negative. Asymptotically, it is guaranteed to be between 0 and 1. A small value indicates low dependence of Y on X, and a high value indicates strong dependence. The statistic is the base for the \code{\link{soFun.test}}.
#' @export
#' @param X Matrix of functional covariates (p by n).
#' @param Y Vector of scalar responses (length n).
#' @param eps Error bound for nearest neighbour search: a value of 0.0 implies exact nearest neighbour search.
#' @param na.rm Remove NAs if TRUE
#' @return The measure of dependence of Y on X.
#' @author Daniel Strenger-Galvis
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
    X <- as.matrix(X[,ok, drop=FALSE])
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
