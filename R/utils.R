logsumexp <- function(x)
{
  m <- max(x)
  m + log(sum(exp(x - m)))
}

msg <- function(s)
{
  time <- format(Sys.time(), "%X")
  message(sprintf("%s %s", time, s))
}

mvn_entropy <- function(sigma)
{
  sigma <- as.matrix(sigma)
  p <- nrow(sigma)
  p/2*(1 + log(2*pi)) + 1/2*log(det(sigma))
}

safe_log <- function(x)
{
  is_zero <- x < 1e-40
  lx <- numeric(length(x))
  lx[!is_zero] <- log(x[!is_zero])
  lx[is_zero] <- -1e4
  lx
}
