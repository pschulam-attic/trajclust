#' Compute posterior for a curve.
#'
#' @param X The design matrix for the curve.
#' @param x The measurement times for the curve.
#' @param y A response vector for the curve.
#' @param K The covariance matrix of the responses.
#' @param model A trajclust model.
#'
#' @export
trajclust_inference <- function(X, x, y, K, model) {
  U <- cbind(1, x)
  Ki <- solve(K)

  z <- numeric(model$num_groups)
  bmean <- matrix(0, length(model$bmean), model$num_groups)
  bcov <- solve(solve(model$bcov) + t(U) %*% Ki %*% U)
  bmean_smooth <- solve(model$bcov) %*% model$bmean

  res_mean <- as.numeric(U %*% model$bmean)
  res_cov <- K + U %*% model$bcov %*% t(U)

  for (i in 1:length(z)) {
    yhat <- as.numeric(X %*% model$beta[, i])
    z[i] <- log(model$theta[i])
    z[i] <- z[i] + mvtnorm::dmvnorm(y - yhat, res_mean, res_cov, log=TRUE)
    bmean[, i] <- bcov %*% (t(U)%*%Ki%*%(y - yhat) + bmean_smooth)
  }

  likelihood <- logsumexp(z)
  z <- exp(z - likelihood)
  list(z=z, bmean=bmean, bcov=bcov, likelihood=likelihood)
}

#' Use a trained model to infer groups and compute likelihood.
#'
#' @param curveset A collection of curves to evaluate.
#' @param model A trained trajclust model.
#'
#' @export
trajclust_full_inference <- function(curveset, model) {
  likelihood <- 0
  z <- matrix(NA, curveset$num_curves, model$num_groups)
  p <- length(model$bmean)
  bmean <- array(NA, c(p, model$num_groups, curveset$num_curves))
  i <- 0

  for (curve in curveset$curves) {
    i <- i + 1
    x <- curve$x
    y <- curve$y
    X <- model$basis(x)
    K <- model$covariance(curve$x)
    inf <- trajclust_inference(X, x, y, K, model)

    z[i, ] <- inf$z
    bmean[, , i] <- inf$bmean
    likelihood <- likelihood + inf$likelihood
  }

  list(z=z, bmean=bmean, likelihood=likelihood)
}
