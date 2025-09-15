# psicor.test -------------------------------------------------------------
#' Psi correlation test for categorical variables.
#' @param Y Vector (length n) of categorical values
#' @param X Matrix (n by q)
#' @param eps Error bound for nearest neighbour search: a value of 0.0 implies exact nearest neighbour search.
#' @param na.rm Remove NAs if TRUE
#' @param permutation If TRUE, the p-value is calculated by permutation. If FALSE, the p-value is calculated based on the asymptotic distribution of the test statistic. The default is FALSE.
#' @param R Number of permutations. Ignored if permutation == FALSE.
#' @details The test statistic is a mean square contingency, similar to the usual chi-squared divergence.
#' It is known to be asymptotically distributed as a chi-squared distribution with K*(K-1) degrees of freedom, where K is the number of levels of Y. The p-value is either calculated based on this distribution or on a resampling procedure, according to the value of the argument \code{permutation}.
#' @return A list containing the test statistic and its p-value.
#' @import RANN
#' @export
#' @author Daniel Strenger-Galvis
#' @references Hörmann, S. and Strenger, D. (2025).  Quantifying and testing dependence to categorical variables. \href{https://arxiv.org/abs/2509.10268}{https://arxiv.org/abs/2509.10268}
#' @examples
#' X <- matrix(runif(100*2), ncol=2)
#' Y <- sin(2*pi*rowSums(X)) >= 0 # A highly non-linear relationship
#' psicor.test(Y,X)
psicor.test <- function(Y, X, eps = 2, na.rm = TRUE, permutation = FALSE, R = 1000){
  if (na.rm == TRUE) {
    # NAs are removed here:
    ok <- complete.cases(Y,X)
    Y <- Y[ok]
    X <- as.matrix(X[ok,])
  }
  Y = factor(Y)
  n = length(Y)
  p.hat = as.vector(table(Y))/n
  K = length(p.hat)

  if(K==n){
    stop("None of the values for Y was observed more than once. Note: This might be a hint that Y actually follows a continuous distribution or that your sample size is too small.")
  }
  if(permutation){
    statistic <- psicor(Y, X, eps = eps)
    permutation_statistics <- numeric(R)
    for(i in 1:R){
      Y.permuted <- sample(Y)
      permutation_statistics[i] <- psicor(Y.permuted, X, eps = eps)
    }
    p.value <- mean(permutation_statistics >= statistic)
  }else{
    neighbor.indices = nn2(X, k=2, eps = eps)$nn.idx[, 2]
    yN = Y[neighbor.indices]
    yN = factor(yN, levels = levels(Y))
    q.hat = as.vector(table(yN))/n
    covariance_matrix = calculate_covariance(p.hat, neighbor.indices)
    covariance_matrix = round(covariance_matrix, 15)


    tab = matrix(0, nrow = K-1, ncol = K-1)
    for(k in 1:(K-1)){
      for(l in 1:(K-1)){
        tab[k,l] = sum(((Y==levels(Y)[k])-p.hat[k])*((yN==levels(Y)[l])-q.hat[l]))/sqrt(n)
      }
    }
    deviations = as.vector(tab)

    inv = solve(covariance_matrix, deviations)
    statistic = deviations%*%inv
    p.value = 1-pchisq(statistic, (K-1)*(K-1))
  }

  return(list(statistic = statistic, p.value = p.value))
}
