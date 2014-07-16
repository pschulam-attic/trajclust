#' @export
plot.trajclust <- function(x, ...) {
  model <- x
  xrange <- model$train_info$xrange
  yrange <- model$train_info$yrange
  yextra <- 1/4 * diff(yrange)

  n <- 100
  x <- seq(xrange[1], xrange[2], length=n)
  X <- model$basis(x)
  y <- X %*% model$beta
  group <- as.factor(seq(model$num_groups))
  trajs <- data.frame(x=x, y=as.numeric(y), group=rep(group, each=n))

  p <- ggplot2::ggplot(trajs)
  p <- p + ggplot2::geom_line(ggplot2::aes(x, y, color=group))
  p <- p + ggplot2::facet_wrap(~ group)
  p <- p + ggplot2::ylim(yrange[1] - yextra, yrange[2] + yextra)
  p
}
