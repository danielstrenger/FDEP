library(RANN)
library(rcompanion)

# psicor -------------------------------------------------------------------------
#' Estimate the psi correlation.
#'
#' Psi correlation is a measure of the amount of conditional dependence between a categorical random variable Y and a random variable Z taking values in a separable metric space given a random variable X taking values in a separable metric space, based on an i.i.d. sample of (Y,Z,X). The coefficient is asymptotically guaranteed to be between 0 and 1.
#' @param Y Vector (length n)
#' @param Z Matrix (n by p)
#' @param X Matrix (n by q), default is NULL
#' @param eps Error bound for nearest neighbour search: a value of 0.0 implies exact nearest neighbour search.
#' @param na.rm Remove NAs if TRUE
#' @details The value returned by psicor can be positive or negative. Asymptotically, it is guaranteed to be between 0 and 1. If X==NULL, the finite sample version is non-negative as well. psicor is invariant under permutation of the classes of Y.
#' @return Psi correlation of Y and Z given X. If X == NULL, this is just a measure of the dependence between Y and Z.
#' @import RANN
#' @import rcompanion
#' @export
#' @author Daniel Strenger-Galvis
#' @references Hörmann, S. and Strenger-Galvis, D. (2025). Quantifying and testing dependence to categorical variables. \href{https://arxiv.org/abs/2509.10268}{https://arxiv.org/abs/2509.10268}
#' @examples
#' Z <- matrix(runif(100*2), ncol=2)
#' Y <- sin(2*pi*rowSums(Z)) >= 0 # A highly non-linear relationship
#' psicor(Y, Z)
psicor <- function(Y, Z, X=NULL, eps = 2, na.rm = TRUE){
  if (na.rm == TRUE) {
    # NAs are removed here:
    ok <- complete.cases(Y,Z,X)
    Y <- Y[ok]
    Z <- as.matrix(Z[ok,])
    if(!is.null(X)){
      X <- as.matrix(X[ok,])
    }
  }
  if(is.null(X)){
    Y <- factor(Y)
    neighbor.indices <- nn2(Z, k=2, eps = eps)$nn.idx[, 2]
    yN <- Y[neighbor.indices]
    yN <- factor(yN, levels = levels(Y))

    if(length(unique(yN)) == 1){
      warning("The nearest neighbour of Y only contains one level. Return NA.")
      # If yN only contains one level, the psi correlation is 0.
      return(NA)
    }

    return(cramerV(Y, yN)[[1]])
  } else {
    psi_YXZ <- psicor(Y, cbind(Z, X))
    psi_YX <- psicor(Y, X)

    return((psi_YXZ - psi_YX) / (1 - psi_YX))
  }
}

calculate_covariance <- function(p, neighbor.indices){
  n = length(neighbor.indices)
  K = length(p)

  neighbor_pairs_number = sum(neighbor.indices[neighbor.indices]==1:n)
  w1 = neighbor_pairs_number
  neighbor_counts = as.vector(table(neighbor.indices))
  common_neighbor_count = sum(neighbor_counts*(neighbor_counts-1))
  w2 = common_neighbor_count

  covariance_structure = array(rep(0, ((K-1))^4), c(K-1, K-1, K-1, K-1))
  for(k1 in 1:(K-1)){
    for(l1 in 1:(K-1)){
      for(k2 in 1:(K-1)){
        for(l2 in 1:(K-1)){
          covariance_structure[k1, l1, k2, l2] <- (p[k1]*p[l1]*p[k2]*p[l2]*(n+w1)+
                                                     p[k1]*p[l1]*p[l2]*(-n*(k1==k2)-w1*(l1==k2))+
                                                     p[k1]*p[l1]*p[k2]*(-w1*(k1==l2)-n*(l1==l2))+
                                                     p[k1]*p[l1]*(n*(k1==k2)*(l1==l2)+w1*(k1==l2)*(l1==k2)))/n
        }
      }
    }
  }
  return(matrix(covariance_structure, nrow = (K-1)*(K-1)))
}
