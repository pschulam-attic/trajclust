#' Create an empty trajclust model.
#'
#' @param k The number of clusters.
#' @param p The number of basis functions.
#' @param basis A function to compute design matrices.
#' @param covariance A function to compute covariance matrices.
#'
#' @export
new_trajclust_model <- function(k, p, basis, covariance)
{
  model <- structure(list(), class="trajclust")
  model$num_groups <- k
  model$num_basis  <- p
  model$theta <- numeric(k)
  model$beta <- matrix(NA, p, k)
  model$basis <- basis
  model$covariance <- covariance
  model
}

#' Create an empty set of trajclust sufficient statistics.
#'
#' @param model A trajclust model.
#'
#' @export
new_trajclust_suffstats <- function(model)
{
  k <- model$num_groups
  p <- model$num_basis
  
  ss <- structure(list(), class="trajclust_ss")
  ss$theta_suffstats <- numeric(k)
  ss$beta_eta1 <- array(0, c(p, p, k))
  ss$beta_eta2 <- array(0, c(p, k))

  ss
}

#' Initialize a trajclust model using data.
#'
#' @param curveset A collection of curves.
#' @param model A trajclust model.
#' 
#' @export
init_trajclust_model <- function(curveset, model)
{
  theta <- runif(model$num_groups)
  model$theta <- theta / sum(theta)

  for (i in 1:model$num_groups)
  {
    curve <- sample(curveset, 1)
    X <- model$basis(curve$x)
    y <- curve$y
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
  total_theta <- sum(ss$theta_suffstats)
  model$theta <- ss$theta_suffstats / total_theta

  for (i in 1:model$num_groups)
  {
    eta1 <- ss$beta_eta1[, , i]
    eta2 <- ss$beta_eta2[, i]
    model$beta[, i] <- solve(eta1, eta2)
  }

  model
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

  covariance <- function(x)
  {
    outer(x, x, kernel)
  }

  covariance
}
