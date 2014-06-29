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
