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

#' rounds a numeric vector probabilistically such that its sum (an integer) is preserved
#' 
#' @param vec numeric vector which sums to an integer
#' @return numeric vector of the same length as vec: rounded vector
round_preserve_sum <- function(vec) {
    
    # if vector already integers, do nothing
    # old implementation -- all.equal(vec - round(vec)) -- fails for large values in vec
    if(isTRUE(all.equal(vec - round(vec), numeric(length(vec))))) {
        return(vec)
    }
    
    # if vector doesn't sum to integer, 
    # round down so it does sum to an integer (ensuring positivity)
    sum_vec <- sum(vec)
    rounded_sum <- round(sum_vec)
    fractional_part <- rounded_sum - sum_vec
    if(!isTRUE(all.equal(rounded_sum - sum_vec, 0))) {
        larger_idx <- which(vec > fractional_part)
        # subtract fractional_part from one of those groups, chosen randomly
        subtract_idx <- sample(larger_idx, 1)
        vec[subtract_idx] <- vec[subtract_idx] - fractional_part
    }
    
    n_round_up <- rounded_sum - sum(floor(vec))
    round_up_idx <- sample.int(length(vec), size = n_round_up, replace = FALSE,
                               prob = vec - floor(vec))
    vec_out <- floor(vec)
    vec_out[round_up_idx] <- ceiling(vec[round_up_idx])
    vec_out
}