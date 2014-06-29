logsumexp <- function(x)
{
  m <- max(x)
  m + log(sum(exp(x - m)))
}

msg <- function(s)
{
  time <- format(Sys.time(), "%X")
  cat(sprintf("%s %s\n", time, s))
}

mvn_entropy <- function(sigma)
{
  sigma <- as.matrix(sigma)
  p <- nrow(sigma)
  p/2*(1 + log(2*pi)) + 1/2*log(det(sigma))
}

mult_entropy <- function(theta)
{
  - sum(log(theta))
}
