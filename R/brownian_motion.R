library(stats)

#' Simulate paths of the standard Brownian Motion on [0,1].
#'
#' @param n Number of paths to simulate.
#' @param p Number of discretisation points (exluding 0) for each path.
#' @return A matrix (p+1 by n) of simulated independent paths of the standard Brownian motion on [0,1].
#' @export
#' @examples
#' x <- brownian.motion(20,200)
#' matplot(seq(0,1,1/200),x, type="l", xlab="t", ylab="x(t)")
brownian.motion<-function(n, p){
  x <- matrix(rnorm(n*(p+1)),ncol=n)/sqrt(p)
  x[1,] <- 0
  x <- apply(x,2,cumsum)
  return(x)
}
