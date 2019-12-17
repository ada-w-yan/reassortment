sensitivity_mutation_rate <- function() {
  hash <- get_hash()
  mutation_prob <- c(2e-3, 2e-5)
  reassortment <- c(FALSE, TRUE)
  par_grid <- expand.grid(mutation_prob = mutation_prob,
                          reassortment = reassortment)
  par_grid$sim_name <- paste0("mutation_prob_", 
                              num2str(par_grid$mutation_prob),
                              "_",
                              reassortment)
  job <- Map(function(x, y, z) obj$enqueue(run_default_pars(reassortment = x,
                                                            sim_name = y,
                                                            hash = hash,
                                                            mutation_prob = z),
                                           par_grid$reassortment,
                                           par_grid$sim_name,
                                           par_grid$mutation_prob))
  job
}

sensitivity_fitness_WM <- function() {
  hash <- get_hash()
  fitness_WM <- c(1, 1.1, 1.5, 1.75, 2)
  reassortment <- c(FALSE, TRUE)
  par_grid <- expand.grid(fitness_WM = fitness_WM,
                          reassortment = reassortment)
  par_grid$sim_name <- paste0("fitness_WM_", 
                              num2str(par_grid$fitness_WM),
                              "_",
                              reassortment)
  job <- Map(function(x, y, z) obj$enqueue(run_default_pars(reassortment = x,
                                                            sim_name = y,
                                                            hash = hash,
                                                            fitness_WM = z),
                                           par_grid$reassortment,
                                           par_grid$sim_name,
                                           par_grid$fitness_WM))
  job
}

sensitivity_fitness_MW <- function() {
  hash <- get_hash()
  reassortment <- c(FALSE, TRUE)
  sim_name <- paste0("fitness_MW_1_", reassortment)
  job <- Map(function(x, y) obj$enqueue(run_default_pars(reassortment = x,
                                                            sim_name = y,
                                                            hash = hash,
                                                            fitness_MW = 1),
                                           reassortment,
                                           sim_name))
  job
}

sensitivity_fitness_MM <- function() {
  hash <- get_hash()
  fitness_MM <- c(.8, .9, 1.1, 1.2)
  reassortment <- c(FALSE, TRUE)
  par_grid <- expand.grid(fitness_MM = fitness_MM,
                          reassortment = reassortment)
  par_grid$sim_name <- paste0("fitness_MM_", 
                              num2str(par_grid$fitness_MM),
                              "_",
                              reassortment)
  job <- Map(function(x, y, z) obj$enqueue(run_default_pars(reassortment = x,
                                                            sim_name = y,
                                                            hash = hash,
                                                            fitness_MM = z),
                                           par_grid$reassortment,
                                           par_grid$sim_name,
                                           par_grid$fitness_MM))
  job
}