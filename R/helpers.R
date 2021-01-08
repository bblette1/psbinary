# Contrast function
h <- function(x, y, contrast) {
  func <- log(x) - log(y)
  if (contrast=="Difference") { func <- x - y }
  if (contrast=="VE") { func <- 1 - x/y }
  return(func)
}

# Define helper function for Imbens-Manski interval computation
f_cstar <- function(c, low, up, maxsig) {
  pnorm(c + (up - low) / maxsig) - pnorm(-c) - 0.95
}
