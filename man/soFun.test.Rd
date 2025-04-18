% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/soFun.R
\name{soFun.test}
\alias{soFun.test}
\title{Scalar on Function test for dependence.}
\usage{
soFun.test(Y, X, eps = 2, na.rm = TRUE, permutation = FALSE, R = 1000)
}
\arguments{
\item{Y}{Vector (length n)}

\item{X}{Matrix (p by n)}

\item{eps}{Error bound for nearest neighbour search: a value of 0.0 implies exact nearest neighbour search.}

\item{na.rm}{Remove NAs if TRUE}

\item{permutation}{If TRUE, the p-value is calculated by permutation. If FALSE, the p-value is calculated based on the asymptotic distribution of the test statistic. The default is FALSE.}

\item{R}{Number of permutations. Ignored if permutation == FALSE.}
}
\value{
A list containing the test statistic and its p-value.
}
\description{
The test is based on the asymptotic distribution of the soFun statistic.
}
\details{
The test statistic is the \code{\link{soFun}} statistic, with an appropriate scaling applied if permuatation == FALSE.
The test is based either on the known asymptotic normal distribution of this statistic or on a resampling procedure, according to the value of the argument \code{permutation}.
}
\examples{
n <- 100
p <- 200
x <- brownian.motion(n,p)
y <- sin(2*pi* apply(x,2,mean))
if(soFun.test(y,x)$p.value < 0.05){
  print("The hypothesis of independence can be rejected at a significance level of 0.05.")
} else{
  print("The hypothesis of independence cannot be rejected at a significance level of 0.05.")
}
}
\references{
Hörmann, S. and Strenger, D. (2025). Azadkia–Chatterjee’s dependence coefficient for infinite dimensional data.
}
\author{
Daniel Strenger
}
