#' Run EM to fit the trajclust model.
#'
#' @param curveset A collection of observed curves.
#' @param model A trajclust model.
#' @param tol Convergence tolerance.
#'
#' @export
run_em <- function(curveset, model, tol=1e-4)
{
  iter <- 0
  convergence <- 1
  likelihood_old <- 0
  
  while (convergence > tol)
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

    convergence <- (likelihood_old - likelihood) / likelihood_old
    likelihood_old <- likelihood

    msg(sprintf("likelihood=%.2f, convergence=%.4f",
                likelihood, convergence))
  }

  model
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
  X <- model$basis(curve$x)
  y <- curve$y
  K <- model$covariance(curve$x)
  inf <- trajclust_inference(X, y, K, model)
  likelihood <- inf$likelihood

  A <- t(X) %*% solve(K)
  eta1 <- A %*% X
  eta2 <- A %*% y

  for (i in 1:model$num_groups)
  {
    z_i <- inf$z[i]
    ss$theta_suffstats[i] <- ss$theta_suffstats[i] + z_i
    ss$beta_eta1[, , i] <- ss$beta_eta1[, , i] + z_i * eta1
    ss$beta_eta2[, i] <- ss$beta_eta2[, i] + z_i * eta2
  }

  list(ss=ss, likelihood=likelihood)
}
