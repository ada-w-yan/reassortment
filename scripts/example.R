run_example <- function() {
  set.seed(2)
  iv <- c(95,0,0,5) #mt,mt  mt,wt  wt,mt   wt,wt
  fitness <- c(1,0,1.1,1)
  n_cells <- 1e4
  pop_size <- 1000
  mutation_prob <- 0 #2e-5: Russell et al. (2012)
  burst_size <- 10
  generations <- 20
  
  results <- simulate_evolution(iv, fitness, burst_size, n_cells, pop_size,
                                generations, mutation_prob,
                                coinfection = TRUE,
                                MOI_dependent_burst_size = TRUE,
                                choose_strain_by_fitness = FALSE,
                                one_strain_produced = FALSE,
                                reassort = TRUE)
  plot_strains(results)
  plot_segments(results)
}
