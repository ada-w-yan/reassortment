#' vapply returning numeric vector of length 1
#' 
#' @param X see vapply
#' @param FUN see vapply
#' @param ... additional arguments to pass to vapply
#' @return output of vapply
#' @export
vnapply <- function(X, FUN, ...) {
  vapply(X, FUN, numeric(1), ...)
}

#' normalises a vector so that it sums to 1
#' 
#' @param vec numeric vector to be normalised
#' @return normalised vector of same length as vec
normalise <- function(vec) {
  
  if(all(vec == 0)) {
    return(vec)
  } else if(any(vec < 0)) {
    stop("vector to be normalised must have all non-negative elements")
  }
  
  vec / sum(vec)
  
}

#' round a (vector of) numbers up with probability equal to their non-integer parts
#' 
#' @param x numeric vector 
#' @return numeric vector of the same length as x
probabilistic_round <- function(x) {
  floor(x) + (x - floor(x) > runif(length(x)))
}