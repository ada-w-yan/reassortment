run_analyses <- function(sim_name, pop_size = 1e6, hash) {
  set.seed(2)
  iv <- c(95,0,0,5) #mt,mt  mt,wt  wt,mt   wt,wt
  fitness <- c(1,0,1.1,1)
  MOI <- .1
  n_cells <- round(pop_size / MOI)
  mutation_prob <- .09#2e-5 #Russell et al. (2012)
  burst_size <- 10
  generations <- 20
  
  sim_wrapper <- function(iv, fitness, burst_size, n_cells, pop_size,
                          generations) {
    function(mutation_prob,
             coinfection,
             MOI_dependent_burst_size,
             choose_strain_by_fitness,
             one_strain_produced,
             reassort, filename, hash)  {
      results <- simulate_evolution(iv, fitness, burst_size, n_cells, pop_size,
                                    generations, mutation_prob,
                                    coinfection,
                                    MOI_dependent_burst_size,
                                    choose_strain_by_fitness,
                                    one_strain_produced,
                                    reassort)
      dir_name <- make_results_folder(filename, hash = hash)
      saveRDS(results, paste0(dir_name, "results.rds"))
      g <- plot_strains(results)
      ggsave(paste0(dir_name, "strains.pdf"), g, width = 10, height = 10, units = "cm")
      g <- plot_segments(results)
      ggsave(paste0(dir_name, "segments.pdf"), g, width = 10, height = 10, units = "cm")
      invisible(results)
    }
  }
  
  sim_fn <- sim_wrapper(iv, fitness, burst_size, n_cells, pop_size,
                        generations)
  
  switch(sim_name,
         "reassort_MOI_dependent" = sim_fn(mutation_prob = 0,
                                                     coinfection = TRUE,
                                                     MOI_dependent_burst_size = TRUE,
                                                     choose_strain_by_fitness = FALSE,
                                                     one_strain_produced = FALSE,
                                                     reassort = TRUE,
                                                     filename = sim_name),
         "reassort_MOI_independent" =   sim_fn(mutation_prob = 0,
         coinfection = TRUE,
         MOI_dependent_burst_size = FALSE,
         choose_strain_by_fitness = FALSE,
         one_strain_produced = FALSE,
         reassort = TRUE,
         filename = sim_name),
  "mutate_MOI_dependent" = sim_fn(mutation_prob,
         coinfection = TRUE,
         MOI_dependent_burst_size = TRUE,
         choose_strain_by_fitness = FALSE,
         one_strain_produced = FALSE,
         reassort = FALSE,
         filename = sim_name),
  "mutate_MOI_independent" = sim_fn(mutation_prob,
         coinfection = TRUE,
         MOI_dependent_burst_size = FALSE,
         choose_strain_by_fitness = FALSE,
         one_strain_produced = FALSE,
         reassort = FALSE,
         filename = sim_name),
  "mutate_MOI_dependent_by_fitness" = sim_fn(mutation_prob,
         coinfection = TRUE,
         MOI_dependent_burst_size = TRUE,
         choose_strain_by_fitness = TRUE,
         one_strain_produced = FALSE,
         reassort = FALSE,
         filename = sim_name),
  "mutate_MOI_independent_by_fitness" = sim_fn(mutation_prob,
         coinfection = TRUE,
         MOI_dependent_burst_size = FALSE,
         choose_strain_by_fitness = TRUE,
         one_strain_produced = FALSE,
         reassort = FALSE,
         filename = sim_name),
  "mutate_MOI_dependent_one_strain" = sim_fn(mutation_prob,
         coinfection = TRUE,
         MOI_dependent_burst_size = TRUE,
         choose_strain_by_fitness = FALSE,
         one_strain_produced = TRUE,
         reassort = FALSE,
         filename = sim_name),
  "mutate_MOI_independent_one_strain" = sim_fn(mutation_prob,
         coinfection = TRUE,
         MOI_dependent_burst_size = FALSE,
         choose_strain_by_fitness = FALSE,
         one_strain_produced = TRUE,
         reassort = FALSE,
         filename = sim_name),
  "mutate_no_coinfection" = sim_fn(mutation_prob,
         coinfection = FALSE,
         MOI_dependent_burst_size = FALSE,
         choose_strain_by_fitness = FALSE,
         one_strain_produced = FALSE,
         reassort = FALSE,
         filename = sim_name))
}

# obj <- setup_cluster()
# analysis_names <- c("reassort_MOI_dependent",
#                     "reassort_MOI_independent",
#                     "mutate_MOI_dependent",
#                     "mutate_MOI_independent",
#                     "mutate_MOI_dependent_by_fitness",
#                     "mutate_MOI_independent_by_fitness",
#                     "mutate_MOI_dependent_one_strain",
#                     "mutate_MOI_independent_one_strain",
#                     "mutate_no_coinfection")
# hash <- get_hash()
# job <- obj$lapply(analysis_names, function(x) run_analyses (x, 1e4, hash))
# semireligious_kob