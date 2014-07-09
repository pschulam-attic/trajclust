#' Create an empty trajclust model.
#'
#' @param G The number of clusters.
#' @param P The number of basis functions.
#' @param basis A function to compute design matrices.
#' @param covariance A function to compute covariance matrices.
#' @param bmean Individual offset prior mean.
#' @param bcov Individual offset prior cov.
#'
#' @export
new_trajclust_model <- function(G, P, basis, covariance, bmean=NULL, bcov=NULL)
{
  if (is.null(bmean))
    bmean <- rep(0, 2)

  if (is.null(bcov))
    bcov <- diag(1e-10, 2)

  model <- structure(list(), class="trajclust")
  model$num_groups <- G
  model$num_basis  <- P
  model$theta <- numeric(G)
  model$beta <- matrix(NA, P, G)
  model$basis <- basis
  model$covariance <- covariance
  model$bmean <- bmean
  model$bcov <- bcov
  model
}

#' Create an empty set of trajclust sufficient statistics.
#'
#' @param model A trajclust model.
#'
#' @export
new_trajclust_suffstats <- function(model)
{
  G <- model$num_groups
  P <- model$num_basis

  ss <- structure(list(), class="trajclust_ss")
  ss$theta_suffstats <- numeric(G)
  ss$beta_eta1 <- array(0, c(P, P, G))
  ss$beta_eta2 <- array(0, c(P, G))

  ss
}

#' Initialize a trajclust model using data.
#'
#' @param curveset A collection of curves.
#' @param model A trajclust model.
#'
#' @export
init_trajclust_model <- function(curveset, model, seed)
{
  set.seed(seed)

  theta <- runif(model$num_groups)
  model$theta <- theta / sum(theta)

  init_groups <- sample(model$num_groups, curveset$num_curves, TRUE)

  for (i in 1:model$num_groups)
  {
    curves <- curveset$curves[init_groups == i]
    x <- do.call("c", lapply(curves, "[[", "x"))
    y <- do.call("c", lapply(curves, "[[", "y"))
    X <- model$basis(x)
    A <- crossprod(X) + diag(1e-2, ncol(X))
    b <- t(X) %*% y
    model$beta[, i] <- solve(A, b)
  }

  model
}

#' Compute maximum likelihood estimate of trajclust model.
#'
#' @param model A trajclust model.
#' @param ss A trajclust sufficient statistics container.
#'
#' @export
trajclust_mle <- function(model, ss)
{
  alpha <- 1     # Pseudo-counts for group probabilities.
  fudge <- 1e-2  # Diagonal value to prevent singular matrices.

  total_theta <- sum(ss$theta_suffstats) + alpha * model$num_groups
  model$theta <- (ss$theta_suffstats + alpha) / total_theta

  for (i in 1:model$num_groups)
  {
    eta1 <- ss$beta_eta1[, , i] + diag(fudge, model$num_basis)
    eta2 <- ss$beta_eta2[, i]
    model$beta[, i] <- solve(eta1, eta2)
  }

  model
}

#' Create a new polynomial basis.
#'
#' @param degree Degree of the polynomial.
#'
#' @export
polynomial_basis <- function(degree)
{
  basis <- function(x)
  {
    powers <- seq(0, degree)
    X <- matrix(0, length(x), length(powers))
    for (i in 1:length(powers))
        X[, i] <- x^powers[i]

    X
  }

  basis
}

#' Create a new B-spline basis.
#'
#' @param xrange The domain of the B-spline model.
#' @param nbasis The number of basis functions to use.
#' @param intercept Should the basis include an intercept?
#' @param degree The degree of the polynomials.
#'
#' @export
bspline_basis <- function(xrange, nbasis, intercept, degree=3)
{
  from <- xrange[1]
  to   <- xrange[2]

  if (intercept)
      nknots <- nbasis - degree + 1
  else
      nknots <- nbasis - degree + 2

  knots <- seq(from, to, length=nknots)
  boundary <- knots[ c(1, nknots)]
  interior <- knots[-c(1, nknots)]

  basis <- function(x)
  {
    X <- splines::bs(x, degree=degree, knots=interior,
                     Boundary.knots=boundary, intercept=intercept)

    unname(X)
  }

  basis
}

#' Create a new diagonal covariance function.
#'
#' @param noise The square root of the measurement variance.
#'
#' @export
diagonal_covariance <- function(noise)
{
  covariance <- function(x1, x2)
  {
    if (missing(x2))
      diag(noise^2, length(x1))
    else
      matrix(0, length(x1), length(x2))
  }

  covariance
}

#' Create a new squared exponential covariance function.
#'
#' @param amp The amplitude of the kernel.
#' @param bw The bandwidth of the kernel.
#' @param noise The amount of observation noise.
#'
#' @export
squared_exp_covariance <- function(amp, bw, noise)
{
  kernel <- function(x, y)
  {
    d <- abs(x - y)
    w <- amp^2 * exp(-1/2 * d^2 / bw^2)
    w <- w + ifelse(x == y, noise^2, 0)
  }

  covariance <- function(x1, x2)
  {
    if (missing(x2))
      x2 <- x1

    outer(x1, x2, kernel)
  }

  covariance
}
