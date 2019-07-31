paper_panels <- function(MOI, fitness_MW = 0, fitness_WM = 1.25, mutation_prob, reassort, pop_size = 1e6, hash, plot_results = FALSE) {
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
  n_replicates <- 100
  run_parallel <- TRUE
  reassort <- as.logical(reassort)
      
  sim_name <- paste(num2str(c(MOI, fitness_MW, fitness_WM, mutation_prob)), collapse = "_")
    sim_name <-  paste0(sim_name, "_", reassort)
      dir_name <- ifelse(missing(hash),
                         make_results_folder(sim_name),
                         make_results_folder(sim_name, hash = hash))
      inputs <- ls()
      inputs <- list_vars_from_environment(inputs)
      saveRDS(inputs, paste0(dir_name, "inputs.rds"))

      sim_func <- function(run_no)  {
       result <- simulate_evolution(iv, fitness, burst_size, n_cells, pop_size,
                                                generations, mutation_prob,
                                                coinfection,
                                                MOI_dependent_burst_size,
                                                choose_strain_by_fitness,
                                                one_strain_produced,
                                                reassort)
         result <-   cbind(result, matrix(run_no, nrow = generations, ncol = 1, dimnames = list(NULL, "run")))
         result
      }
      results <- parLapply_wrapper(run_parallel, seq_len(n_replicates), sim_func)
        results <-   do.call(rbind, results)

      saveRDS(results, paste0(dir_name, "results.rds"))
      if(plot_results) {
        g <- plot_multirun_strains(results)
        ggsave(paste0(dir_name, "strains.pdf"), g, width = 10, height = 10, units = "cm")
        saveRDS(g, "strains.rds")
        g <- plot_multirun_segments(results)
        ggsave(paste0(dir_name, "segments.pdf"), g, width = 10, height = 10, units = "cm")
        saveRDS(g, "segments.rds")
      }
      invisible(results)
}
if(FALSE) {
paper_panels_separate_seed <- function(MOI, fitness_MW = 0, fitness_WM = 1.25, mutation_prob, reassort, pop_size = 1e6, hash, seed, plot_results = TRUE) {
  set.seed(seed)
  iv <- c(95,0,0,5) #mt,mt  mt,wt  wt,mt   wt,wt
  fitness <- c(1, fitness_MW, fitness_WM, 1)
  n_cells <- round(pop_size / MOI)
  burst_size <- 10
  generations <- 20
  coinfection <- TRUE
  MOI_dependent_burst_size <- TRUE
  choose_strain_by_fitness <- FALSE
  one_strain_produced <- FALSE
  n_replicates <- 12
  run_parallel <- TRUE
  reassort <- as.logical(reassort)
      
  sim_name <- paste(num2str(c(MOI, fitness_MW, fitness_WM, mutation_prob, seed)), collapse = "_") %>%
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

      saveRDS(results, paste0(dir_name, "results.rds"))
      g <- plot_multirun_strains(results)
      ggsave(paste0(dir_name, "strains.pdf"), g, width = 10, height = 10, units = "cm")
      saveRDS(g, "strains.rds")
      g <- plot_multirun_segments(results)
      ggsave(paste0(dir_name, "segments.pdf"), g, width = 10, height = 10, units = "cm")
      saveRDS(g, "segments.rds")
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
  reassort <- TRUE
  
  sim_name <- paste(num2str(MOI), mutation_prob, "constant_cells", sep = "_")
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
  
  saveRDS(results, paste0(dir_name, "results.rds"))
  g <- plot_multirun_strains(results)
  ggsave(paste0(dir_name, "strains.pdf"), g, width = 10, height = 10, units = "cm")
  saveRDS(g, "strains.rds")
  g <- plot_multirun_segments(results)
  ggsave(paste0(dir_name, "segments.pdf"), g, width = 10, height = 10, units = "cm")
  saveRDS(g, "segments.rds")
  invisible(results)
}

get_results <- function(MOI, fitness_MW = 0, fitness_WM = 1.25, mutation_prob, reassort, pop_size = 1e6, hash, job_idx) {
  
  sim_name <- paste(num2str(as.numeric(c(MOI, fitness_MW, fitness_WM, mutation_prob))), collapse = "_") %>%
    paste0("_", trimws(reassort))
  dir_name <- ifelse(missing(hash),
                     make_results_folder(sim_name),
                     make_results_folder(sim_name, hash = hash))
  
  job <- obj$task_get(jobs$ids[as.numeric(job_idx)])
  if(job$status() == "COMPLETE") {
    results <- job$result()
    saveRDS(results, paste0(dir_name, "results.rds"))
    g <- plot_multirun_strains(results)
    ggsave(paste0(dir_name, "strains.pdf"), g, width = 10, height = 10, units = "cm")
    saveRDS(g, "strains.rds")
    g <- plot_multirun_segments(results)
    ggsave(paste0(dir_name, "segments.pdf"), g, width = 10, height = 10, units = "cm")
    saveRDS(g, "segments.rds")
  }
  invisible(NULL)
}

plot_end_change_MOI <- function() {
  MOI <- 10^seq(-3, 1)
  fitness_MW <- 0
  fitness_WM <- 1.25
  mutation_prob <- 2e-4
  reassort <- TRUE
  # using results where the number of virions rather than the number of cells is held constant
  extract_prop_MM_by_MOI <- function(MOI) {
    filename <- paste(num2str(c(MOI, fitness_MW, fitness_WM, mutation_prob)), collapse = "_") %>%
    paste0("_", reassort, "/results.rds") %>%
      paste0("results/", .)
    results <- readRDS(filename)
    n_runs <- max(results[,"run"])
    n_gen <- nrow(results) / n_runs
    results <- results[seq_len(nrow(results)) %% n_gen == 0, "MM"] %>%
      quantile(prob = c(0.025, 0.975))
    results
  }
  if(FALSE) {
  prop_MM_by_MOI <- t(vapply(MOI, extract_prop_MM_by_MOI, double(2))) %>%
    as.data.frame
  prop_MM_by_MOI$MOI <- MOI
  saveRDS(prop_MM_by_MOI, "results/prop_MM_by_MOI.rds")
  }
  
  prop_MM_by_MOI <- readRDS("results/prop_MM_by_MOI.rds")
  
  make_label <- function(x) TeX(paste0("$10^{", x, "}$"))
  x_labels <- sapply(log10(MOI), make_label)
  g <- ggplot(prop_MM_by_MOI, aes(x = MOI, ymin = `2.5%`, ymax = `97.5%`)) +
    geom_errorbar() +
    scale_x_log10(breaks = MOI, labels = x_labels) +
    theme_bw() +
    coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    ylab("Proportion of MM") +
    theme(text = element_text(size = 20))
  ggsave("results/prop_MM_by_MOI.pdf", g, width = 10, height = 10, units = "cm")
  saveRDS(g, "prop_MM_by_MOI_plot.rds")
  invisible(prop_MM_by_MOI)
}

}
if(FALSE) {
  for (i in 4:10) {
    a <- paper_panels_separate_seed(MOI = 1, fitness_MW = 0, fitness_WM = 1.25, mutation_prob = 2e-4, reassort = TRUE, pop_size = 1e6, hash = "57f2fe3", seed = i)
  }
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

  # get_results
  pars$hash <- "2909c83"
  pars$job_idx <- seq_len(nrow(pars))
  apply_named_args(pars, 1, get_results)
  
  
  obj <- setup_cluster(n_cores = 5)
  job <- obj$enqueue_bulk(pars, paper_panels)
  # contributable_paintedladybutterfly (4 to go)
  
  diff_MOI_idx <- seq(4, 6)
  pars <- rbind(pars_reassort_only[diff_MOI_idx,],
                pars_reassort_mutate[diff_MOI_idx,])
  pars$n_cells <- 1e6
  pars$hash <- get_hash()
  pars <- pars[,c("MOI", "mutation_prob", "n_cells", "hash")]
  job <- obj$enqueue_bulk(pars, paper_panels_constant_n_cells)
  # intracranial_lightningbug (done)
  
  pars_reassort_only <- data.frame(MOI = 1e-1,
                                   fitness_MW = 0,
                                   fitness_WM = fitness_WM,
                                   reassort = TRUE,
                                   mutation_prob = 0)
  
  pars_reassort_mutate <- pars_reassort_only
  pars_reassort_mutate$mutation_prob <- mutation_prob
  pars <- rbind(pars_reassort_only,
                pars_reassort_mutate)
  
  pars$pop_size <- 1e6
  pars$hash <- get_hash()
  
  obj <- setup_cluster(n_cores = 5)
  job <- obj$enqueue_bulk(pars, paper_panels)
  # fel_motmot (done)
  
  colnames(pars)[6] <- "n_cells"
  pars <- pars[,c("MOI", "mutation_prob", "n_cells", "hash")]
  job <- obj$enqueue_bulk(pars, paper_panels_constant_n_cells)
  # debatable_hapuka (failed)
  


}