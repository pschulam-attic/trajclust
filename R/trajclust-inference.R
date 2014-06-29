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
