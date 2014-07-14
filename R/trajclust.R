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
#' @param ninit The number of random restarts to use.
#' @param seed A random seed for reproducibility.
#' @param verbose Logical flog indicating whether to print convergence
#' information.
#'
#' @export
trajclust <- function(x, y, id, ngroups, xrange=range(x), nbasis,
                      amp, bw, noise, ninit=5, bmean=NULL, bcov=NULL,
                      seed=1, verbose=TRUE)
{
  set.seed(seed)
  init_seeds <- sample(100, ninit, replace=FALSE)

  curveset <- make_curveset(x, y, id)
  basis <- bspline_basis(xrange, nbasis, TRUE)
  covariance <- squared_exp_covariance(amp, bw, noise)

  model <- new_trajclust_model(ngroups, nbasis, basis, covariance, bmean, bcov)
  model$train_info <- list()
  model$train_info$xrange <- curveset$xrange
  model$train_info$yrange <- curveset$yrange

  fit_models <- vector("list", ninit)
  likelihoods <- numeric(ninit)

  for (i in 1:ninit)
  {
    if (verbose) cat(banner(i), "\n")

    s <- init_seeds[i]
    m <- init_trajclust_model(curveset, model, seed=s)
    em <- run_em(curveset, m, verbose=verbose)
    m <- em$model
    m$train_info$likelihood <- em$likelihood
    fit_models[[i]] <- m
  }

  best <- which.max(likelihoods)
  fit_models[[best]]
}

banner <- function(iter)
{
  head <- paste0(rep("#", 76), collapse="")
  head <- paste0("\n", head)
  text <- paste("###", sprintf("Starting trajclust restart %d\n", iter))
  paste(head, text, sep="\n")
}
