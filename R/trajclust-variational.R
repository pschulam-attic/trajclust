#' Run variational EM to fit the trajclust model. The model must have
#' been initialized with \code{bmean} and \code{bcov}.
#'
#' This procedure is used to estimate the model when individual curves
#' can have linear offsets.
#'
#' @param curveset A collection of observed curves.
#' @param model A trajclust model.
#' @param tol Convergence tolerance.
#'
#' @export
run_var_em <- function(curveset, model, tol=1e-8, maxiter=1e4, verbose=TRUE)
{
  # Initialize iteration variables.

  iter <- 0
  convergence <- 1
  likelihood_old <- 0

  while (convergence > tol & iter < maxiter)
  {
    # Initialize iteration.

    iter <- iter + 1
    likelihood <- 0
    ss <- new_trajclust_suffstats(model)

    # Collect sufficient statistics and likelihood for each curve.

    for (curve in curveset$curves)
    {
      estep <- curve_var_e_step(curve, model, ss)
      ss <- estep$ss
      likelihood <- likelihood + estep$likelihood
    }

    # Use sufficient statistics to compute trajclust MLE.

    model <- trajclust_mle(model, ss)

    # Update convergence criterion.

    convergence <- abs((likelihood_old - likelihood) / likelihood_old)
    likelihood_old <- likelihood

    if (verbose)
        msg(sprintf("iter=%04d, likelihood=%.2f, convergence=%.8f",
                    iter, likelihood, convergence))
  }

  list(model=model, likelihood=likelihood)
}

#' Compute sufficient statistics and likelihood using a variational
#' approximation.
#'
#' @param curve The curve to process.
#' @param model A trajclust model.
#' @param ss A trajclust sufficient statistics container.
#'
#' @export
curve_var_e_step <- function(curve, model, ss)
{
  # Extract curve observations.

  X <- model$basis(curve$x)
  x <- curve$x
  y <- curve$y
  K <- model$covariance(curve$x)

  # Compute posterior over latent variables.

  inf <- trajclust_var_inference(X, x, y, K, model)
  likelihood <- inf$likelihood

  # Update sufficient statistics.

  z <- inf$z
  bmean <- inf$bmean

  A <- t(X) %*% solve(K)
  eta1 <- A %*% X
  eta2 <- A %*% (y - cbind(1, x) %*% bmean)

  for (i in 1:model$num_groups)
  {
    ss$theta_suffstats[i] <- ss$theta_suffstats[i] + z[i]
    ss$beta_eta1[, , i] <- ss$beta_eta1[, , i] + z[i] * eta1
    ss$beta_eta2[, i] <- ss$beta_eta2[, i] + z[i] * eta2
  }

  list(ss=ss, likelihood=likelihood, z=z, bmean=bmean)
}

#' Compute posterior group and offset for a curve.
#'
#' @param X The design matrix for the curve.
#' @param x The measurement times for the curve.
#' @param y The observed measurements for the curve.
#' @param K The measurement covariance matrix.
#' @param model A trajclust model.
#' @param tol Convergence tolerance for variational inference.
#'
#' @export
trajclust_var_inference <- function(X, x, y, K, model, tol=1e-8)
{
  iter <- 0
  convergence <- 1
  likelihood_old <- 0

  n <- length(x)
  Xb <- cbind(1, x)
  Ki <- solve(K)
  yb <- y - X %*% model$beta
  precision <- t(Xb) %*% Ki %*% Xb
  projection <- t(Xb) %*% Ki %*% yb

  # Initialize variational distributions; var_bcov is constant.

  var_z <- numeric(model$num_groups)
  var_bmean <- numeric(2)
  var_bcov <- model$bcov

  while (convergence > tol)
  {
    iter <- iter + 1

    # Update var_z.

    ebsq <- var_bcov + var_bmean %*% t(var_bmean)

    for (i in 1:model$num_groups)
    {
      zi <- safe_log(model$theta[i])
      zi <- zi - n/2 * log(2*pi)
      zi <- zi - 1/2 * log(det(K))
      zi <- zi - 1/2 * yb[, i] %*% Ki %*% yb[, i]
      zi <- zi - 1/2 * sum(diag(t(Xb) %*% Ki %*% Xb %*% ebsq))
      zi <- zi + projection[, i] %*% var_bmean
      var_z[i] <- zi
    }
    var_z <- exp(var_z - logsumexp(var_z))

    # Update var_bmean.

    var_bmean <- solve(model$bcov) %*% model$bmean

    for (i in 1:model$num_groups)
    {
      zi <- var_z[i]
      var_bmean <- var_bmean + zi*projection[, i]
    }

    var_bmean <- solve(solve(var_bcov) + precision, var_bmean)

    likelihood <- trajclust_elbo(X, x, y, K, var_z, var_bmean, var_bcov, model)

    convergence <- abs((likelihood_old - likelihood) / likelihood_old)
    likelihood_old <- likelihood

    ## msg(sprintf("elbo=%.2f, convergence=%.8f", likelihood, convergence))
  }

  list(z=var_z, bmean=var_bmean, bcov=var_bcov, likelihood=likelihood)
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
  likelihood <- likelihood + t(model$bmean) %*% bcovi %*% var_bmean
  likelihood <- likelihood + mvn_entropy(var_bcov)

  # E[prior(z)] + H(var(z))
  likelihood <- likelihood + sum(var_z * safe_log(model$theta))
  likelihood <- likelihood - sum(var_z * safe_log(var_z))

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

#' Use a trained model to infer groups and offsets and compute
#' likelihood.
#'
#' @param curveset A collection of curves to evaluate.
#' @param model A trained trajclust model.
#'
#' @export
trajclust_full_var_inference <- function(curveset, model)
{
  likelihood <- 0
  z <- matrix(NA, curveset$num_curves, model$num_groups)
  bmean <- matrix(NA, curveset$num_curves, 2)
  i <- 0

  for (curve in curveset$curves)
  {
    i <- i + 1
    X <- model$basis(curve$x)
    x <- curve$x
    y <- curve$y
    K <- model$covariance(curve$x)
    inf <- trajclust_var_inference(X, x, y, K, model)

    z[i, ] <- inf$z
    bmean[i, ] <- inf$bmean
    likelihood <- likelihood + inf$likelihood
  }

  list(z=z, bmean=bmean, likelihood=likelihood)
}
