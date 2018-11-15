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
  # old implementation probably safer, but slower...
  # if(isTRUE(all.equal(vec - round(vec), numeric(length(vec))))) {
  #     return(vec)
  # }
  # ever older implementation -- all.equal(vec, round(vec)) -- fails for large values in vec
  if(all(vec - round(vec) - numeric(length(vec)) == 0)) {
    return(vec)
  }
  
  # if vector doesn't sum to integer, 
  # round down so it does sum to an integer (ensuring positivity)
  sum_vec <- sum(vec)
  rounded_sum <- round(sum_vec)
  fractional_part <- rounded_sum - sum_vec
  if(any(rounded_sum - sum_vec != 0)) {
    # old implementation probably safer, but slower...
    # if(!isTRUE(all.equal(rounded_sum - sum_vec, 0))) {
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

#' get short hash of current commit
#' 
#' @param repo_dir repo directory
#' @return character vector containing short hash of current commit
#' @export
get_hash <- function(repo_dir) {
  if(missing(repo_dir)) {
    repo_dir <- getwd()
  }
  command <- paste("git -C", repo_dir, "rev-parse --short HEAD")
  hash <- system(command, intern = TRUE)
  hash
}

#' make folder to store results and (optional) store hash of repo
#' 
#' @param folder_name folder name
#' @param base_dir directory in which to make folder
#' @param hash string -- make file with this name in folder
#' @return name of created directory
#' @importFrom magrittr %>% %T>%
#' @export
make_results_folder <- function(folder_name, base_dir, hash) {
  if(missing(base_dir)) {
    base_dir <- paste0(file.path(getwd(), "results"), "/")
  }
  
  # make folder (warn if exists)
  dir_name <- paste0(base_dir, folder_name, "/") %T>%
    dir.create(., recursive = TRUE)
  # make file named after current has of repo to which base_dir belongs
  if(!missing(hash)) {
    system(paste0("touch ", dir_name, hash))
  }
  
  dir_name
}

#' gather variables with given names from an environment into a list
#' 
#' @param var_names character vector of variable names to gather
#' @param envir environment in which to find parameters
#' @return list of variables with names var_names
#' @export
list_vars_from_environment <- function(var_names, envir = parent.frame()) {
  env_vars <- as.list(envir)
  env_vars <- get_vars_from_list_with_check(env_vars, var_names)
  env_vars
}

#' extract variables from list, throwing an error if they are not found
#' 
#' @param x list of variables
#' @param var_names character vector containing names of variables to extract
#' @return list of selected variables
#' @export
get_vars_from_list_with_check <- function(x, var_names) {
  missing_vars <- var_names[!(var_names %in% names (x))]
  if(length(missing_vars) > 0) {
    stop(paste0("variables missing from list: ", paste0(missing_vars, collapse = " ")))
  }
  x <- x[var_names]
  x
}

#' turn number(s) into string for filename, replacing decimal points with "point"
#' 
#' @param x numeric vector
#' @return character vector of same length as x
#' @export
num2str <- function(x) {
  gsub("\\.", "point", as.character(x))
}

#' wrapper for parLapply for cluster
#'
#' @param run_parallel: logical: if TRUE, use parLapply, else use lapply
#' @param x first argument of lapply
#' @param fun second argument of lapply
#' @return output arguments of lapply
#' @export
parLapply_wrapper <- function(run_parallel,x,fun,...){
  if(run_parallel){
    sys_info <- Sys.info()
    if(sys_info[[1]] == "Windows"){
      parallel::parLapply(cl = NULL, x, fun, ...)
    } else {
      parallel::mclapply(x, fun, ..., mc.cores = length(x))
    }
  } else {
    lapply(x, fun, ...)
  }
}