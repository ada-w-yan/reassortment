change_fitness_MOI_mutation_prob <- function(fitness_WM, MOI, mutation_prob, reassort, pop_size = 1e6, hash) {
  set.seed(2)
  iv <- c(95,0,0,5) #mt,mt  mt,wt  wt,mt   wt,wt
  fitness <- c(1, 0, fitness_WM, 1)
  n_cells <- round(pop_size / MOI)
  burst_size <- 10
  generations <- 20
  coinfection <- TRUE
  MOI_dependent_burst_size <- TRUE
  choose_strain_by_fitness <- FALSE
  one_strain_produced <- FALSE
  n_replicates <- 10
  run_parallel <- TRUE
  reassort <- as.logical(reassort)
      
  sim_name <- paste(num2str(c(fitness_WM, MOI, mutation_prob)), collapse = "_")
      dir_name <- ifelse(missing(hash),
                         make_results_folder(sim_name),
                         make_results_folder(sim_name, hash = hash))
      inputs <- ls()
      inputs <- list_vars_from_environment(inputs)
      saveRDS(inputs, paste0(dir_name, "inputs.rds"))

      sim_func <- function(run_no)  {
       simulate_evolution(iv, fitness, burst_size, n_cells, pop_size,
                                                generations, mutation_prob,
                                                coinfection,
                                                MOI_dependent_burst_size,
                                                choose_strain_by_fitness,
                                                one_strain_produced,
                                                reassort) %>%
          cbind(., matrix(run_no, nrow = generations, ncol = 1, dimnames = list(NULL, "run")))
      }
      results <- parLapply_wrapper(run_parallel, seq_len(n_replicates), sim_func) %>%
        do.call(rbind, .)

      # results <- readRDS(paste0(dir_name, "results.rds"))
      g <- plot_multirun_strains(results)
      ggsave(paste0(dir_name, "strains.pdf"), g, width = 10, height = 10, units = "cm")
      g <- plot_multirun_segments(results)
      ggsave(paste0(dir_name, "segments.pdf"), g, width = 10, height = 10, units = "cm")
      invisible(results)
}

if(FALSE) {
  fitness_WM <- c(1.1, 2, 5, 10)
  MOI <- c(0.1, 1, 5)
  mutation_prob <- 10^seq(-5, -3)
  reassort_pars <- expand.grid(fitness_WM = fitness_WM, MOI = MOI)
  reassort_pars$mutation_prob <- 0
  reassort_pars$reassort <- TRUE
  reassort_pars$pop_size <- 1e6
  # apply_named_args(reassort_pars, 1, change_fitness_MOI_mutation_prob)
  # reassort_pars$hash <- get_hash()
  # change_fitness_MOI_mutation_prob(fitness_WM, MOI[1], 0, TRUE, 1e3)
  # obj <- setup_cluster(n_cores = 1)
  # job <- obj$enqueue_bulk(reassort_pars, change_fitness_MOI_mutation_prob)
  
  mutation_pars <- expand.grid(fitness_WM = fitness_WM, MOI = MOI, mutation_prob = mutation_prob)
  mutation_pars$reassort <- FALSE
  mutation_pars$pop_size <- 1e6
  mutation_pars$hash <- get_hash()
  obj <- setup_cluster(n_cores = 5)
  job <- obj$enqueue_bulk(mutation_pars, change_fitness_MOI_mutation_prob)
}