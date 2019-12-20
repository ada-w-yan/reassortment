default_wrapper <- function() {
  hash <- get_hash()
  reassort <- c(FALSE, TRUE)
  par_grid <- expand.grid(reassort = reassort)
  par_grid$sim_name <- paste0("default_", 
                              par_grid$reassort)
  job <- obj$enqueue_bulk(par_grid, run_default_pars)
  job
}

sensitivity_mutation_rate <- function() {
  hash <- get_hash()
  mutation_prob <- c(2e-3, 2e-5)
  reassort <- c(FALSE, TRUE)
  par_grid <- expand.grid(mutation_prob = mutation_prob,
                          reassort = reassort)
  par_grid$sim_name <- paste0("mutation_prob_", 
                              num2str(par_grid$mutation_prob),
                              "_",
                              par_grid$reassort)
  job <- obj$enqueue_bulk(par_grid, run_default_pars)
  job
}

sensitivity_fitness_WM <- function() {
  hash <- get_hash()
  fitness_WM <- c(1, 1.1, 1.5, 1.75, 2)
  reassort <- c(FALSE, TRUE)
  par_grid <- expand.grid(fitness_WM = fitness_WM,
                          reassort = reassort)
  par_grid$sim_name <- paste0("fitness_WM_", 
                              num2str(par_grid$fitness_WM),
                              "_",
                              par_grid$reassort)
  job <- obj$enqueue_bulk(par_grid, run_default_pars)
  job
}

sensitivity_fitness_MW <- function() {
  hash <- get_hash()
  reassort <- c(FALSE, TRUE)
  par_grid <- expand.grid(fitness_MW = 1,
                          reassort = reassort)
  par_grid$sim_name <- paste0("fitness_MW_1_",
                              par_grid$reassort)
  job <- obj$enqueue_bulk(par_grid, run_default_pars)
  job
}

sensitivity_fitness_MM <- function() {
  hash <- get_hash()
  fitness_MM <- c(.8, .9, 1.1, 1.2)
  reassort <- c(FALSE, TRUE)
  par_grid <- expand.grid(fitness_MM = fitness_MM,
                          reassort = reassort)
  par_grid$sim_name <- paste0("fitness_MM_", 
                              num2str(par_grid$fitness_MM),
                              "_",
                              par_grid$reassort)
  job <- obj$enqueue_bulk(par_grid, run_default_pars)
  job
}

sensitivity_MOI <- function() {
  hash <- get_hash()
  MOI <- c(.1, 10)
  reassort <- c(FALSE, TRUE)
  par_grid <- expand.grid(MOI = MOI,
                          reassort = reassort)
  par_grid$sim_name <- paste0("MOI_", 
                              num2str(par_grid$MOI),
                              "_",
                              par_grid$reassort)
  job <- obj$enqueue_bulk(par_grid, run_default_pars)
  job
}

sensitivity_burst_size <- function() {
  hash <- get_hash()
  burst_size <- c(1, 100)
  reassort <- c(FALSE, TRUE)
  par_grid <- expand.grid(burst_size = burst_size,
                          reassort = reassort)
  par_grid$sim_name <- paste0("MOI_", 
                              num2str(par_grid$burst_size),
                              "_",
                              par_grid$reassort)
  job <- obj$enqueue_bulk(par_grid, run_default_pars)
  job
}

sensitivity_MOI_dependent_burst_size <- function() {
  hash <- get_hash()
  par_grid <- expand.grid(MOI_dependent_burst_size = FALSE,
                          reassort = reassort)
  par_grid$sim_name <- paste0("MOI_dependent_burst_size_FALSE_", 
                              par_grid$reassort)
  job <- obj$enqueue_bulk(par_grid, run_default_pars)
  job
}