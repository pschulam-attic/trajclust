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
#' @param verbose Logical flog indicating whether to print convergence
#' information.
#'
#' @export
trajclust <- function(x, y, id, ngroups, xrange=range(x), nbasis,
                      amp, bw, noise, ninit=5, bmean=NULL, bcov=NULL,
                      verbose=TRUE)
{
  curveset <- make_curveset(x, y, id)
  basis <- bspline_basis(xrange, nbasis, TRUE)
  covariance <- squared_exp_covariance(amp, bw, noise)

  model <- new_trajclust_model(ngroups, nbasis, basis, covariance, bmean, bcov)
  model$train_info <- list()
  model$train_info$xrange <- curveset$xrange
  model$train_info$yrange <- curveset$yrange

  model <- init_trajclust_model(curveset, model)
  em <- run_em(curveset, m, verbose=verbose)
  model <- em$model
  model$train_info$likelihood <- em$likelihood

  model
}

banner <- function(iter)
{
  head <- paste0(rep("#", 76), collapse="")
  head <- paste0("\n", head)
  text <- paste("###", sprintf("Starting trajclust restart %d\n", iter))
  paste(head, text, sep="\n")
}
