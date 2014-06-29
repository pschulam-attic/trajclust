#' Compute posterior for a curve.
#'
#' @param X A design matrix for the curve.
#' @param y A response vector for the curve.
#' @param K The covariance matrix of the responses.
#' @param model A trajclust model.
#'
#' @export
trajclust_inference <- function(X, y, K, model)
{
  z <- numeric(model$num_groups)

  for (i in 1:length(z))
  {
    yhat <- X %*% model$beta[, i]
    z[i] <- mvtnorm::dmvnorm(y, yhat, K, log=TRUE)
    z[i] <- z[i] + log(model$theta[i])
  }

  likelihood <- logsumexp(z)
  z <- exp(z - likelihood)
  list(z=z, likelihood=likelihood)
}

#' Use a trained model to infer groups and compute likelihood.
#'
#' @param curveset A collection of curves to evaluate.
#' @param model A trained trajclust model.
#'
#' @export
trajclust_full_inference <- function(curveset, model)
{
  likelihood <- 0
  z <- matrix(NA, curveset$num_curves, model$num_groups)
  i <- 0

  for (curve in curveset$curves)
  {
    i <- i + 1
    X <- model$basis(curve$x)
    y <- curve$y
    K <- model$covariance(curve$x)
    inf <- trajclust_inference(X, y, K, model)

    z[i, ] <- inf$z
    likelihood <- likelihood + inf$likelihood
  }

  list(z=z, likelihood=likelihood)
}
