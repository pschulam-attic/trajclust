#' Run EM to fit the trajclust model.
#'
#' @param curveset A collection of observed curves.
#' @param model A trajclust model.
#' @param tol Convergence tolerance.
#'
#' @export
run_em <- function(curveset, model, tol=1e-8, maxiter=1e4, verbose=TRUE)
{
  iter <- 0
  convergence <- 1
  likelihood_old <- 0

  while (convergence > tol & iter < maxiter)
  {
    iter <- iter + 1

    likelihood <- 0
    ss <- new_trajclust_suffstats(model)

    for (curve in curveset$curves)
    {
      estep <- curve_e_step(curve, model, ss)
      ss <- estep$ss
      likelihood <- likelihood + estep$likelihood
    }

    model <- trajclust_mle(model, ss)
    ## capture.output(print(model), file=sprintf("trajclust-%03d.txt", iter))

    convergence <- abs((likelihood - likelihood_old) / likelihood_old)
    likelihood_old <- likelihood

    if (verbose)
      msg(sprintf("likelihood=%.2f, convergence=%.8f",
                  likelihood, convergence))
  }

  list(model=model, likelihood=likelihood)
}

#' Compute sufficient statistics and likelihood.
#'
#' @param curve The curve to process.
#' @param model A trajclust model.
#' @param ss A trajclust sufficient statistics container.
#'
#' @export
curve_e_step <- function(curve, model, ss)
{
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

  for (i in 1:model$num_groups)
  {
    ss$theta_suffstats[i] <- ss$theta_suffstats[i] + z[i]
    ss$beta_eta1[, , i] <- ss$beta_eta1[, , i] + z[i] * eta1
    eta2 <- A %*% (y - U %*% bmean[, i])
    ss$beta_eta2[, i] <- ss$beta_eta2[, i] + z[i] * eta2
    ss$beta_cov_ss[, , i] <- ss$beta_cov_ss[, , i] + z[i] * beta_cov_ss
  }

  list(ss=ss, likelihood=likelihood)
}
