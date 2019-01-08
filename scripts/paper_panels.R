library(magrittr)
paper_panels <- function(MOI, fitness_MW = 0, fitness_WM = 1.25, mutation_prob, reassort, pop_size = 1e6, hash) {
  set.seed(2)
  iv <- c(95,0,0,5) #mt,mt  mt,wt  wt,mt   wt,wt
  fitness <- c(1, fitness_MW, fitness_WM, 1)
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
      
  sim_name <- paste(num2str(c(MOI, fitness_MW, mutation_prob)), collapse = "_") %>%
    paste0("_", reassort)
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

paper_panels_constant_n_cells <- function(MOI, mutation_prob, n_cells = 1e6, hash) {
  set.seed(2)
  iv <- c(95,0,0,5) #mt,mt  mt,wt  wt,mt   wt,wt
  fitness <- c(1, 0, 1.25, 1)
  pop_size <- round(n_cells * MOI)
  burst_size <- 10
  generations <- 20
  coinfection <- TRUE
  MOI_dependent_burst_size <- TRUE
  choose_strain_by_fitness <- FALSE
  one_strain_produced <- FALSE
  n_replicates <- 10
  run_parallel <- TRUE
  reassort <- as.logical(reassort)
  
  sim_name <- paste(num2str(MOI), reassort, "constant_cells", collapse = "_")
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
  mutation_prob <- 2e-4
  fitness_WM <- 1.25
  pars_reassort_only <- data.frame(MOI = c(1, 1, 1, 1e-3, 1e-2, 10),
                     fitness_MW = c(0, 1, 0, 0, 0, 0),
                     fitness_WM = fitness_WM - c(0, 0, fitness_WM - 1, 0, 0, 0))
  pars_reassort_only$reassort <- TRUE
  pars_reassort_only$mutation_prob <- 0

  pars_reassort_mutate <- pars_reassort_only
  pars_reassort_mutate$mutation_prob <- mutation_prob
  
  pars_mutate_only <- pars_reassort_mutate[1,]
  pars_mutate_only$reassort <- FALSE
  
  pars <- rbind(pars_reassort_only,
                pars_reassort_mutate,
                pars_mutate_only)
  
  pars$pop_size <- 1e6
  pars$hash <- get_hash()

  obj <- setup_cluster(n_cores = 5)
  job <- obj$enqueue_bulk(pars, paper_panels)
  
  diff_MOI_idx <- seq(4, 6)
  pars <- rbind(pars_reassort_only[diff_MOI_idx,],
                pars_reassort_mutate[diff_MOI_idx,])
  pars$n_cells <- 1e6
  pars$hash <- get_hash()
  job <- obj$enqueue_bulk(pars, paper_panels_constant_n_cells)
}