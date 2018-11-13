run_analyses <- function() {
  set.seed(2)
  iv <- c(95,0,0,5) #mt,mt  mt,wt  wt,mt   wt,wt
  fitness <- c(1,0,1.1,1)
  n_cells <- 1e4
  pop_size <- 1000
  mutation_prob <- 2e-5 #Russell et al. (2012)
  burst_size <- 10
  generations <- 20
  
  sim_wrapper <- function(iv, fitness, burst_size, n_cells, pop_size,
                          generations) {
    function(mutation_prob,
             coinfection,
             MOI_dependent_burst_size,
             choose_strain_by_fitness,
             one_strain_produced,
             reassort, filename)  {
      results <- simulate_evolution(iv, fitness, burst_size, n_cells, pop_size,
                                    generations, mutation_prob,
                                    coinfection,
                                    MOI_dependent_burst_size,
                                    choose_strain_by_fitness,
                                    one_strain_produced,
                                    reassort)
      saveRDS(results, paste0("results/", filename, ".rds"))
      g <- plot_strains(results)
      ggsave(paste0("results/", filename, "_strains.pdf"), g, width = 10, height = 10, units = "cm")
      g <- plot_segments(results)
      ggsave(paste0("results/", filename, "_segments.pdf"), g, width = 10, height = 10, units = "cm")
      invisible(results)
    }
  }
  
  sim_fn <- sim_wrapper(iv, fitness, burst_size, n_cells, pop_size,
                        generations)
  
  sim_fn(mutation_prob = 0,
         coinfection = TRUE,
         MOI_dependent_burst_size = TRUE,
         choose_strain_by_fitness = FALSE,
         one_strain_produced = FALSE,
         reassort = TRUE,
         filename = "results_reassort_MOI_dependent")
  
  # sim_fn(mutation_prob = 0,
  #        coinfection = TRUE,
  #        MOI_dependent_burst_size = FALSE,
  #        choose_strain_by_fitness = FALSE,
  #        one_strain_produced = FALSE,
  #        reassort = TRUE,
  #        filename = "results_reassort_MOI_independent")
  # 
  # sim_fn(mutation_prob,
  #        coinfection = TRUE,
  #        MOI_dependent_burst_size = TRUE,
  #        choose_strain_by_fitness = FALSE,
  #        one_strain_produced = FALSE,
  #        reassort = FALSE,
  #        filename = "results_mutate_MOI_dependent")
  # 
  # sim_fn(mutation_prob,
  #        coinfection = TRUE,
  #        MOI_dependent_burst_size = FALSE,
  #        choose_strain_by_fitness = FALSE,
  #        one_strain_produced = FALSE,
  #        reassort = FALSE,
  #        filename = "results_mutate_MOI_independent")
  # 
  # sim_fn(mutation_prob,
  #        coinfection = TRUE,
  #        MOI_dependent_burst_size = TRUE,
  #        choose_strain_by_fitness = TRUE,
  #        one_strain_produced = FALSE,
  #        reassort = FALSE,
  #        filename = "results_mutate_MOI_dependent_by_fitness")
  # 
  # sim_fn(mutation_prob,
  #        coinfection = TRUE,
  #        MOI_dependent_burst_size = FALSE,
  #        choose_strain_by_fitness = TRUE,
  #        one_strain_produced = FALSE,
  #        reassort = FALSE,
  #        filename = "results_mutate_MOI_independent_by_fitness")
  # 
  # sim_fn(mutation_prob,
  #        coinfection = TRUE,
  #        MOI_dependent_burst_size = TRUE,
  #        choose_strain_by_fitness = FALSE,
  #        one_strain_produced = TRUE,
  #        reassort = FALSE,
  #        filename = "results_mutate_MOI_dependent_one_strain")
  # 
  # sim_fn(mutation_prob,
  #        coinfection = TRUE,
  #        MOI_dependent_burst_size = FALSE,
  #        choose_strain_by_fitness = FALSE,
  #        one_strain_produced = TRUE,
  #        reassort = FALSE,
  #        filename = "results_mutate_MOI_independent_one_strain")
  # 
  # sim_fn(mutation_prob,
  #        coinfection = FALSE,
  #        MOI_dependent_burst_size = FALSE,
  #        choose_strain_by_fitness = FALSE,
  #        one_strain_produced = FALSE,
  #        reassort = FALSE,
  #        filename = "results_mutate_no_coinfection")
}