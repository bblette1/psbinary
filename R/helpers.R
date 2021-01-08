# Contrast function
h <- function(x, y, contrast) {
  func <- log(x) - log(y)
  if (contrast=="Difference") { func <- x - y }
  if (contrast=="VE") { func <- 1 - x/y }
  return(func)
}
