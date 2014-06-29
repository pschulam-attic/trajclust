#' @export
plot.trajclust <- function(model, x)
{
  n <- length(x)
  X <- model$basis(x)
  y <- X %*% model$beta
  group <- as.factor(seq(model$num_groups))
  trajs <- data.frame(x=x, y=as.numeric(y), group=rep(group, each=n))

  p <- ggplot2::ggplot(trajs)
  p <- p + ggplot2::geom_line(aes(x, y, color=group))
  p <- p + facet_wrap(~ group)
  p
}
