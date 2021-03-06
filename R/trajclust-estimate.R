#' Run EM to fit the trajclust model.
#'
#' @param curveset A collection of observed curves.
#' @param model A trajclust model.
#' @param tol Convergence tolerance.
#' @param maxiter The maximum number of EM iterations.
#' @param verbose Logical flag indicating whether to print convergence
#' information.
#'
#' @export
run_em <- function(curveset, model, tol=1e-8, maxiter=1e3, verbose=TRUE) {
  iter <- 0
  convergence <- 1
  likelihood_old <- 0

  while (convergence > tol & iter < maxiter) {
    iter <- iter + 1

    likelihood <- 0
    ss <- new_trajclust_suffstats(model)

    for (curve in curveset$curves) {
      estep <- curve_e_step(curve, model, ss)
      ss <- estep$ss
      likelihood <- likelihood + estep$likelihood
    }

    model <- trajclust_mle(model, ss)

    convergence <- abs((likelihood - likelihood_old) / likelihood_old)
    likelihood_old <- likelihood

    if (verbose)
      msg(sprintf("likelihood=%.2f, convergence=%.8f", likelihood, convergence))
  }

  nparam <- length(model$theta) + length(model$beta) + length(model$bmean) + length(model$bcov) + 3
  aic <- -2 * likelihood + 2 * nparam
  aic_c <- aic + 2 * nparam * (nparam + 1) / (curveset$num_points - nparam - 1)
  bic <- -2 * likelihood + nparam * log(curveset$num_points)

  list(model=model, likelihood=likelihood, aic=aic, aic_c=aic_c, bic=bic)
}

#' Compute sufficient statistics and likelihood.
#'
#' @param curve The curve to process.
#' @param model A trajclust model.
#' @param ss A trajclust sufficient statistics container.
#'
#' @export
curve_e_step <- function(curve, model, ss) {
  x <- curve$x
  y <- curve$y
  X <- model$basis(curve$x)
  U <- cbind(1, x)
  K <- model$covariance(curve$x)
  inf <- trajclust_inference(X, x, y, K, model)
  likelihood <- inf$likelihood

  z <- inf$z
  bmean <- inf$bmean

  A <- t(X) %*% solve(K)
  eta1 <- A %*% X
  eta2 <- A %*% y
  beta_cov_ss <- t(X) %*% solve(K + U%*%model$bcov%*%t(U)) %*% X
  ss$bcov_suffstats <- ss$bcov_suffstats + inf$bcov

  for (i in 1:model$num_groups) {
    ss$theta_suffstats[i] <- ss$theta_suffstats[i] + z[i]
    ss$beta_eta1[, , i] <- ss$beta_eta1[, , i] + z[i] * eta1
    eta2 <- A %*% (y - U %*% bmean[, i])
    ss$beta_eta2[, i] <- ss$beta_eta2[, i] + z[i] * eta2
    ss$beta_cov_suffstats[, , i] <- ss$beta_cov_suffstats[, , i] + z[i] * beta_cov_ss
    ss$bcov_suffstats <- ss$bcov_suffstats + z[i] * outer(bmean[, i], bmean[, i])
  }

  list(ss=ss, likelihood=likelihood)
}
