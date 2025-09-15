library(doSNOW)
library(snow)
library(pbmcapply)
utils::globalVariables(c("i")) # to avoid code checker complaining about i having no binding

# psicor.vs -------------------------------------------------------------------------
#' Variable selection based on psi correlation.
#'
#' A variable selection algorithm on psi correlation \link{psicor}.
#' @param Y Vector (length n)
#' @param X Matrix (n by p)
#' @param threshold Threshold for psi correlation: the algorithm stops if the conditional psi correlation is below
#' this value.
#' @param num_variables Maximum number of features to select. If the number of features selected
#' reaches this value, the algorithm stops.
#' @param cores Number of cores to use for parallel computation. If cores > 1, the function uses parallel computation.
#' @param verbose If TRUE, the function prints progress messages
#' @details A forward stepwise variable selection algorithm using on psi correlation \link{psicor} at each step.
#' It terminates when the conditional psi correlation of any possible additional
#' variable and the response given the already selected variables is below \code{threshold}
#' or if \code{num_variables} variables are selected. To reduce the computation time,
#' in each step the conditional psi correlation is calculated for all remaining variables in parallel if \code{cores}>1.
#' @return Psi correlation of Y and Z given X. If X == NULL, this is just a measure of the dependence between Y and Z.
#' @import RANN
#' @import rcompanion
#' @import doSNOW
#' @import pbmcapply
#' @import snow
#' @import foreach
#' @export
#' @author Daniel Strenger-Galvis
#' @references @references Hörmann, S. and Strenger-Galvis, D. (2025). Quantifying and testing dependence to categorical variables. \href{https://arxiv.org/abs/2509.10268}{https://arxiv.org/abs/2509.10268}
#' @examples
#' library(kernlab)
#' data(spam)
#' X <- spam[, -58]
#' Y <- spam[, 58]
#' psicor.vs(Y, X, threshold = 0, cores = parallel::detectCores())
psicor.vs <- function(Y,X, threshold = 0, num_variables = Inf, cores = 1, verbose = TRUE){
  if(cores > 1){
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
    on.exit(stopCluster(cl))
  }

  if(verbose){
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
  }else{
    opts <- list()
  }


  psi <- Inf
  psi_unconditional_old <- 0
  selected <- c()
  while((psi > threshold) && (length(selected) < ncol(X))){
    if(length(selected) == 0){
      excluded <- 1:ncol(X)
    }else{
      excluded <- (1:ncol(X))[-selected]
    }
    if(verbose){
      print(paste("Selecting variable", length(selected) + 1))
    }
    pb <- progressBar(min = 1, max = length(excluded), style="ETA")
    coefficients <- foreach(i = excluded, .combine = c, .options.snow = opts, .export = c("psicor"),
                            .packages = c("rcompanion", "RANN")) %dopar% {
      if(length(selected) == 0){
        X_i <- matrix(X[, i], ncol = 1)
        psicor(Y, X_i)
      }else{
        psicor(Y, X[, c(selected, i)])
      }
    }
    i_max <- which.max(coefficients)
    psi_unconditional <- coefficients[i_max]

    psi <- (psi_unconditional - psi_unconditional_old) / (1 - psi_unconditional_old)

    if(psi > threshold){
      selected <- c(selected, excluded[i_max])
      psi_unconditional_old <- psi_unconditional
    }

    if(length(selected)== num_variables){
      if(verbose){
        print(paste("Stopping because num_variables =", num_variables, "has been reached."))
      }
      break
    }
  }

  return(list(selected = selected,
              psi = psi_unconditional))
}
