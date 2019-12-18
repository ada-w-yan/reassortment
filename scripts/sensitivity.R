sensitivity_mutation_rate <- function() {
  hash <- get_hash()
  mutation_prob <- c(2e-3, 2e-5)
  reassort <- c(FALSE, TRUE)
  par_grid <- expand.grid(mutation_prob = mutation_prob,
                          reassort = reassort)
  par_grid$sim_name <- paste0("mutation_prob_", 
                              num2str(par_grid$mutation_prob),
                              "_",
                              reassort)
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
                              reassort)
  job <- obj$enqueue_bulk(par_grid, run_default_pars)
  job
}

sensitivity_fitness_MW <- function() {
  hash <- get_hash()
  reassort <- c(FALSE, TRUE)
  par_grid <- expand.grid(fitness_MW = 1,
                          reassort = reassort)
  par_grid$sim_name <- paste0("fitness_MW_", 
                              num2str(par_grid$fitness_MW_),
                              "_",
                              reassort)
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
                              reassort)
  job <- obj$enqueue_bulk(par_grid, run_default_pars)
  job
}