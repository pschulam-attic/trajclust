#' Create curves from "narrow" data.
#'
#' @param x A vector of measurement times.
#' @param y A vector of measurements corresponding to \code{x}.
#' @param curve_id Curve IDs for measurements in \code{x} and
#' \code{y}.
#'
#' @export
make_curveset <- function(x, y, curve_id)
{
  curve_data <- data.frame(x=x, y=y, id=curve_id)
  curveset <- structure(list(), class="curveset")
  curveset$curves <- by(curve_data, curve_data$id, make_curve)
  curveset$num_curves <- length(curveset$curves)
  curveset$num_points <-
      sum(vapply(curveset$curves, "[[", integer(1), "num_points"))
  curveset
}


#' Create an individual curve from a data.frame.
#'
#' @param curve_data A \code{data.frame} containing measurement times
#' \code{x}, measurements \code{y}, and a curve ID \code{id}.
#'
#' @export
make_curve <- function(curve_data)
{
  curve <- structure(list(), class="curve")
  curve$id <- curve_data$id[1]
  curve$x  <- curve_data$x
  curve$y  <- curve_data$y
  curve$num_points <- nrow(curve_data)
  curve
}
