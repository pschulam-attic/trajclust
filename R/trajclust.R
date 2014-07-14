#' Cluster trajectories.
#'
#' @param x A vector of measurement times.
#' @param y A vector of measurement values.
#' @param id A vector of trajectory IDs for each (x, y) pair.
#' @param ngroups The number of clusters to estimate.
#' @param xrange The lower and upper bounds of the measurement times.
#' @param nbasis The number of basis functions to use.
#' @param amp The amplitude of the correlation model.
#' @param bw The bandwidth of the correlation model.
#' @param noise The observation standard deviation.
#' @param bmean Normal mean parameter for parametric random effects.
#' @param bcov Normal covariance parameter for parametric random effects.
#' @param verbose Logical flog indicating whether to print convergence
#' information.
#'
#' @export
trajclust <- function(x, y, id, ngroups, xrange=range(x), nbasis,
                      amp, bw, noise, bmean=NULL, bcov=NULL, verbose=TRUE)
{
  curveset <- make_curveset(x, y, id)
  basis <- bspline_basis(xrange, nbasis, TRUE)

  hyper <- expand.grid(amp=amp, bw=bw, noise=noise)
  iter_msg <- paste0("(%0", nchar(as.character(nrow(hyper))), "d)")
  models <- NULL

  for (i in 1:nrow(hyper)) {
    a <- hyper$amp[i]
    b <- hyper$bw[i]
    n <- hyper$noise[i]
    iter_str <- sprintf(iter_msg, i)
    msg(sprintf("%s Fitting with hyperparameters (amp=%.01f, bw=%.01f, noise=%.01f).", iter_str, a, b, n))
    covariance <- squared_exp_covariance(a, b, n)
    model <- new_trajclust_model(ngroups, nbasis, basis, covariance, bmean, bcov)
    model$train_info$xrange <- curveset$xrange
    model$train_info$yrange <- curveset$yrange
    model$train_info$amp <- a
    model$train_info$bw <- b
    model$train_info$noise <- n
    model <- init_trajclust_model(curveset, model)
    em <- run_em(curveset, model, tol=1e-1, verbose=FALSE)
    model <- em$model
    model$train_info$likelihood <- em$likelihood
    models <- c(models, list(model))
  }

  cat("\n")
  msg("Tuning model with highest likelihood.")
  likelihoods <- vapply(models, function(m) m$train_info$likelihood, numeric(1))
  model <- models[[which.max(likelihoods)]]

  em <- run_em(curveset, model, verbose=verbose)
  em$model
}
