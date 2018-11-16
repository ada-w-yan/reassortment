#' Cluster setup
#'
#' This function sets up an object linking to the DIDE cluster in a very crude way.
#' @param n_cores number of cluster cores
#' @return a didehpc::queue_didehpc object
#' @export
setup_cluster <- function(n_cores = 1){
    user_options <- get_user_options()
    setwd(user_options$wd)
    sources <- list.files(path = "scripts", 
                            pattern = "\\.R", 
                            recursive = FALSE,
                            full.names = TRUE)
    do.call(options, user_options$cluster_options)
    src <- provisionr::package_sources(local = user_options$package_dir, expire = 1e-10)
    ## Setup contexts
    context::context_log_start()
    root <- "contexts"
    packages <- list(attached = c("reassortment", "ggplot2", "magrittr", "parallel"),
                     loaded = c("dplyr", "tidyr"))
    ctx <- context::context_save(packages=packages,path=root, sources=sources,package_sources=src)
    config <- didehpc::didehpc_config(cores = n_cores)
    ## Submit setup to cluster
    obj1 <- didehpc::queue_didehpc(ctx, config)
    return(obj1)
}

#' specify user options for cluster
#' 
#' @return a list with the working directory on the network drive; the directory
#' in which the package code sits; and options for didehpc
#' @export
get_user_options <- function() {
  # wd is the working directory to run the cluster job from. 
  # This should be the user's network home drive eg. "~/net/home/reassortment"
    list(wd = "~/net/home/reassortment/",
         package_dir = "~/git_repos/reassortment/",
         cluster_options = list(didehpc.username = "ayan",
                                didehpc.cluster = "fi--didemrchnb"))

}