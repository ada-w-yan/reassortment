#' Simulate evolution of mixture by reassortment
#' 
#' @param iv numeric vector of length n_strains = 4.  
#' Initial proportions of [mt, mt], [mt, wt], [wt, mt], [wt, wt].
#' @param fitness numeric vector of length n_strains = 4.
#' Relative fitness of the four strains.
#' @param burst_size numeric vector of length 1.
#' Burst size from a cell infected with one virion with fitness 1.
#' @param n_cells numeric vector of length 1.  Number of cells.
#' @param pop_size numeric vector of length 1. Initial number of virions.
#' @param max_pop_size numeric vector of length 1. maximum number of virions.
#' @param generations numeric vector of length 1. Number of generations for
#' which to run model.
#' @return matrix of dim c(generations, n_strains).  Proportion of virions of
#' each strain for each generation.
#' @importFrom magrittr %>%
#' @export
run_reassortment_model <- function(iv, fitness, burst_size, n_cells, pop_size, max_pop_size, generations) {
  n_strains <- 4
  # assign strains to initial virus population
  strain <- sample.int(n_strains, 
                       size=population_size, 
                       replace=TRUE, prob=normalise(iv))
  strain_segments <- matrix(c(1, 1, 1, 0, 0, 1, 0, 0), ncol = 2, byrow = TRUE)
  assign_cells <- function(virus_popn) {
    cbind(virus_popn, sample.int(n_cells, nrow(virus_popn), replace=TRUE))
  }
  # decide which virion is in which cell (assume Poisson distribution)
  virus_popn <- strain_segments[strain,] %>%
    assign_cells
  colnames(virus_popn) <- c("mt_pb1", "mt_pa", "cell_no")
  
  segments_to_strain <- function(idx) {
    -(2 * idx[1] + idx[2] + 1) + 5
  }
  segments_to_strain_wrapper <- function(virus_popn) {
    apply(virus_popn[, c("mt_pb1", "mt_pa"), drop = FALSE], 1, segments_to_strain)
  }
  sum_strains <- function(strains) {
    vnapply(seq_len(n_strains), function(x) sum(strains == x))
  }
  
  results <- matrix(0, nrow = generations, ncol = n_strains)
  results[1,] <- (segments_to_strain_wrapper(virus_popn) %>% sum_strains)
  for(generation in 2:generations)  {
    infected_cells <- unique(virus_popn[,"cell_no"])
    make_new_popn <- function(infected_cell) {
      viruses_in_cell <- virus_popn[,"cell_no"]==infected_cell
      strains_in_cell <- segments_to_strain_wrapper(virus_popn[viruses_in_cell,,drop = FALSE])
      # assume actual burst size is a function of average burst size and fitnesses of each virus in the cell
      # and is deterministic
      burst_size_from_cell <- probabilistic_round(burst_size*sum(fitness * sum_strains(strains_in_cell)),0)  
      # assume replication efficiency of each segment is the same, so fitness affects packaging only
      prop_mt_in_cell <- colSums(virus_popn[viruses_in_cell,c("mt_pb1", "mt_pa"), drop = FALSE]) / sum(viruses_in_cell)
      new_popn <- vapply(prop_mt_in_cell, 
                         function(x) rbinom(burst_size_from_cell, 1, x),
                         numeric(burst_size_from_cell))
      new_popn
    }

    virus_popn <- lapply(infected_cells, make_new_popn) %>%
      do.call(rbind, .) %>%
      assign_cells # assumes each new generation infects the same number of cells --
    # probably ok if effective reproduction number = 1
    
    colnames(virus_popn) <- c("mt_pb1", "mt_pa", "cell_no")
    # cap virus population at maximum population
    # if effective reproduction number = 1, virus population is constant -- set
    # initial population equal to max population
    if(nrow(virus_popn) > max_pop_size) {
      virus_popn <- virus_popn[sample.int(nrow(virus_popn), max_pop_size),]
    }
    results[generation,] <- (segments_to_strain_wrapper(virus_popn) %>% sum_strains)
  }
  results <- apply(results, 1, normalise)
  
  segment_to_strain_name <- function(idx) {
    paste0(idx, collapse = "") %>%
      gsub("1", "M", .) %>%
      gsub("0", "W", .)
  }
  
  rownames(results) <- apply(strain_segments, 1, segment_to_strain_name)
  t(results)
}

#' Simulate evolution of mixture by reassortment
#' 
#' @param iv numeric vector of length n_strains = 4.  
#' Initial proportions of [mt, mt], [mt, wt], [wt, mt], [wt, wt].
#' @param fitness numeric vector of length n_strains = 4.
#' Relative fitness of the four strains.
#' @param burst_size numeric vector of length 1.
#' Burst size from a cell infected with one virion with fitness 1.
#' @param pop_size numeric vector of length 1. Initial number of virions.
#' @param generations numeric vector of length 1. Number of generations for
#' which to run model.
#' @param mutation_prob probability of mutation for each new virion.  Assume same
#' for pb1 and pa.
#' @return matrix of dim c(generations, n_strains).  Proportion of virions of
#' each strain for each generation.
#' @importFrom magrittr %>%
#' @export
run_mutation_model_no_coinfection <- function(iv, fitness, burst_size, pop_size,
                                              generations, mutation_prob) {
  n_strains <- 4
  # assign strains to initial virus population
  popn_by_strain <- as.numeric(rmultinom(1, size=population_size, prob=normalise(iv))
  burst_size_by_strain_per_virus <- burst_size * fitness
  results <- matrix(0, nrow = generations, ncol = n_strains)
  results[1,] <- popn_by_strain
  
  mutation_matrix <- matrix(c(0, mutation_prob, mutation_prob, mutation_prob^2,
                            mutation_prob, 0, mutation_prob^2, mutation_prob,
                            mutation_prob, mutation_prob^2, 0, mutation_prob,
                            mutation_prob^2, mutation_prob, mutation_prob, 0),
                            nrow = n_strains, byrow = TRUE)
  diag(mutation_matrix) <- 1 - rowSums(mutation_matrix)
  
  for(generation in 2:generations)  {
    burst_size <- probabilistic_round(burst_size_by_strain_per_virus * popn_by_strain) %>%
      matrix(., ncol = 1) %>% 
      matmult(mutation_matrix, .) %>%
      as.numeric
    # should be selecting from, rather than sampling with probability...
    popn_by_strain <- rmultinom(1, population_size, normalise(burst_size)) %>%
      as.numeric
    results[generation,] <- popn_by_strain
  }
  
  results <- apply(results, 1, normalise)
  
  rownames(results) <- c("MM", "MW", "WM", "WW")
  t(results)
}

#' plot results of run_reassortment_model by strain
#' 
#' @param results output of run_reassortment_model  
#' @import ggplot2
#' @export
plot_strains <- function(results) {
  results <- as.data.frame(results)
  results$gen <- seq_len(nrow(results))
  results_melt <- reshape2::melt(results, id.var = "gen")
  ggplot(results_melt, aes(x = gen, y = value, color = variable, group = variable)) +
    geom_line() +
    theme_bw() +
    coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    xlab("Generation") +
    ylab("Proportion of viral load")
}

#' plot results of run_reassortment_model by segment
#' 
#' @param results output of run_reassortment_model  
#' @import ggplot2
#' @export
plot_segments <- function(results) {
  results <- as.data.frame(results)
  results$gen <- seq_len(nrow(results))
  results$mt_pb1 <- results$MM + results$MW
  results$mt_pa <- results$MM + results$WM
  results_melt <- reshape2::melt(results[,c("gen", "mt_pb1", "mt_pa")], id.var = "gen")
  ggplot(results_melt, aes(x = gen, y = value, color = variable, group = variable)) +
    geom_line() +
    theme_bw() +
    coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    xlab("Generation") +
    ylab("Proportion of viral load")
  
}