library(RANN)
library(stats)


# compute_covariance -------------------------------------------------------------------------
#' Compute the conditional covariance matrix given X.
#'
#' Computes the conditional covariance matrix given X of the empirical deviations from independence under the null hypothesis.
#' Only intended for use in the \link{psicor.test} function.
#' @param p Vector (length K)
#' @param neighbor.indices Vector (length n)
#' @return The conditional covariance matrix.
#' @import RANN
#' @import stats
#' @export
#' @author Daniel Strenger-Galvis
#' @references Hörmann, S. and Strenger-Galvis, D. (2025). Quantifying and testing dependence to categorical variables. \href{https://arxiv.org/abs/2509.10268}{https://arxiv.org/abs/2509.10268}
#' @examples
#' library(RANN)
#' x <- brownian.motion(100, 200)
#' p <- c(0.3, 0.7)
#'
#' neighbor.indices <- nn2(x, k=2)$nn.idx[, 2]
#' sigma <- compute_covariance(p, neighbor.indices)
compute_covariance <- function(p, neighbor.indices){
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
