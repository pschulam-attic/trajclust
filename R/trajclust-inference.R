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

trajclust_var_inference <- function(X, x, y, K, model, tol=1e-8)
{
  iter <- 0
  convergence <- 1
  likelihood_old <- 0

  Xb <- cbind(1, x)
  Ki <- solve(K)
  precision <- t(Xb) %*% Ki %*% Xb
  projection <- t(y - X %*% model$beta) %*% Ki %*% Xb

  # Initialize variational distributions
  var_z <- rep(1, model$num_groups) / model$num_groups
  var_bmean <- model$bmean
  var_bcov <- mdoel$bcov

  while (convergence > tol)
  {
    iter <- iter + 1

    likelihood <- trajclust_elbo(X, x, y, K, var_z, var_bmean,
                                 var_bcov, model)

    # Update z_var

    exp_outer <- var_bcov + var_bmean %*% t(var_bmean)
    
    for (i in 1:model$num_groups)
    {
      zi <- log(model$theta[i])
      zi <- zi - 1/2*sum(diag(precision %*% exp_outer))
      zi <- zi - projection[, i] %*% var_bmean
      z_var[i] <- zi
    }
    z_var <- z_var / sum(z_var)

    # Update bmean var_bmean and var_bcov

    var_bcov <- solve(solve(model$bcov + precision))
    var_bmean <- solve(model$bcov) %*% model$bmean

    for (i in 1:model$num_groups)
    {
      zi <- z_var[i]
      var_bmean <- var_bmean + zi*projection[, i]
    }

    var_bmean <- solve(var_bcov, var_bmean)

    convergence <- (likelihood_old - likelihood) / likelihood_old
    likelihood_old <- likelihood

    msg(sprintf("elbo=%.2f, convergence=%.8f",
                likelihood, convergence))
  }
}

#' Compute the ELBO of the trajclust model with the given variational
#' distributions.
#'
#' @param X A design matrix for a curve.
#' @param x The measurement times for a curve.
#' @param y The observed measurements for a curve.
#' @param K The covariance matrix of the observed measurements.
#' @param var_z The variational distribution over groups.
#' @param var_bmean The variational mean over curve offsets.
#' @param var_bcov The variational covariance over curve offsets.
#' @param model A trajclust model.
#'
#' @export
trajclust_elbo <- function(X, x, y, K, var_z, var_bmean, var_bcov, model)
{
  n <- nrow(X)
  bcovi <- solve(model$bcov)
  exp_outer <- var_bcov + var_bmean %*% t(var_bmean)
  likelihood <- 0

  # E[prior(b)] + H(var(b))
  likelihood <- likelihood - log(2*pi)
  likelihood <- likelihood - 1/2*log(det(model$bcov))
  likelihood <- likelihood - 1/2*t(model$bmean) %*% bcovi %*% model$bmean
  likelihood <- likelihood - 1/2*sum(diag(bcovi %*% exp_outer))
  likelihood <- likelihood - t(model$bmean) %*% bcovi %*% var_bmean
  likelihood <- likelihood + mvn_entropy(var_bcov)

  # E[prior(z)] + H(var(z))
  likelihood <- likelihood + sum(var_z * log(model$theta))
  likelihood <- likelihood + mult_entropy(var_z)

  # E[likelihood(y)]
  likelihood <- likelihood - n/2*log(2*pi)
  likelihood <- likelihood - 1/2*log(det(K))
  Ki <- solve(K)
  Xb <- cbind(1, x)

  for (i in 1:model$num_groups)
  {
    zi <- var_z[i]
    res <- y - X %*% model$beta[, i]

    likelihood <- likelihood - 1/2 * zi * t(res) %*% Ki %*% res
    likelihood <- likelihood - 1/2 * sum(diag(t(Xb) %*% Ki %*% Xb %*% exp_outer))
    likelihood <- likelihood + zi * t(res) %*% Ki %*% Xb %*% var_bmean
  }

  likelihood
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
