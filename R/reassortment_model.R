#' Simulate evolution of mixture by reassortment and/or mutation, for many replicates
#' @param MOI multiplicity of infection = number of virions/number of cells
#' @param iv numeric vector of length n_strains = 4.  
#' Initial proportions of [mt, mt], [mt, wt], [wt, mt], [wt, wt].
#' @param fitness_MW fitness of the PB1 mutant
#' @param fitness_WM fitness of the PA mutant
#' @param fitness_MM fitness of the double mutant
#' @param mutation_prob probability of mutation for each new virion.  Assume same
#' for pb1 and pa.
#' @param reassort logical. If TRUE, random reassortment of segments occurs
#'  before packaging.  If FALSE, each strain in the cell effectively produces
#'  pre-packaged viruses
#' @param pop_size numeric vector of length 1. Number of virions.
#' @param burst_size numeric vector of length 1.
#' Burst size from a cell infected with one virion with fitness 1.
#' @param sim_name character vector for naming of output file
#' @param hash character vector: git hash for tagging run
#' @return matrix of dim c(generations * n_replicates, n_strains + 1).  
#' Concatenated simulation results for n_replicates simulations with the
#' above model parameters.
#' @importFrom magrittr %>%
#' @export
run_default_pars <- function(MOI = 1,
                             iv = c(95,0,0,5),
                             fitness_MW = 0.01,
                             fitness_WM = 1.25,
                             fitness_MM = 1,
                             mutation_prob = 2e-4,
                             reassort = TRUE, 
                             pop_size = 1e6, 
                             burst_size = 10,
                             MOI_dependent_burst_size = TRUE,
                             choose_strain_by_fitness = FALSE,
                             sim_name,
                             hash) {
  set.seed(2)
  fitness <- c(1, fitness_MW, fitness_WM, fitness_MM)
  n_cells <- round(pop_size / MOI)
  generations <- 20
  coinfection <- TRUE
  one_strain_produced <- FALSE
  n_replicates <- 100
  run_parallel <- TRUE
  reassort <- as.logical(reassort)
  
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
  results <- parLapply_wrapper(run_parallel, detectCores(), sim_func)
  results <-   do.call(rbind, results)
  
  saveRDS(results, paste0(dir_name, "results.rds"))
  invisible(results)
}

#' Simulate evolution of mixture by mutation, with no co-infection
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
#' @export
run_mutation_model_no_coinfection <- function(iv, fitness, burst_size, pop_size,
                                              generations, mutation_prob) {
  # equivalent to allowing co-infection, but only one virion is selected to be
  # productive and that virion dictates the burst size... or not?  Because we're
  # saying that the virions are uniformly spread out across cells rather than
  # being Poisson distributed
  n_strains <- 4
  # assign strains to initial virus population
  popn_by_strain <- round_preserve_sum(normalise(iv) * pop_size)
  burst_size_by_strain_per_virus <- burst_size * fitness
  results <- matrix(0, nrow = generations, ncol = n_strains)
  results[1,] <- popn_by_strain
  
  mutation_matrix <- make_mutation_matrix(mutation_prob)
  
  for(generation in 2:generations)  {
    burst_size <- probabilistic_round(burst_size_by_strain_per_virus * popn_by_strain)
    burst_size <- matrix(burst_size, ncol = 1)
    # deterministic mutation
    burst_size <- (mutation_matrix %*% burst_size)
    burst_size <- as.numeric(burst_size)
    # to do: should be selecting from, rather than sampling with probability...
    popn_by_strain <- rmultinom(1, pop_size, normalise(burst_size))
    popn_by_strain <- as.numeric(popn_by_strain)
    results[generation,] <- popn_by_strain
  }
  
  results <- apply(results, 1, normalise)
  
  rownames(results) <- c("MM", "MW", "WM", "WW")
  t(results)
}

#' Simulate evolution of mixture by reassortment and/or mutation
#' 
#' @param iv numeric vector of length n_strains = 4.  
#' Initial proportions of [mt, mt], [mt, wt], [wt, mt], [wt, wt].
#' @param fitness numeric vector of length n_strains = 4.
#' Relative fitness of the four strains.
#' @param burst_size numeric vector of length 1.
#' Burst size from a cell infected with one virion with fitness 1.
#' @param n_cells numeric vector of length 1.  Number of cells.
#' @param pop_size numeric vector of length 1. Number of virions.
#' @param generations numeric vector of length 1. Number of generations for
#' which to run model.
#' @param mutation_prob probability of mutation for each new virion.  Assume same
#' for pb1 and pa.
#' @param coinfection logical.  If TRUE, cells can be infected by more than one
#' virion.  If FALSE, each cell can only be infected by one virion.
#' @param MOI_dependent_burst_size logical.  If TRUE, for a cell infected with a
#' single strain, the burst size scales with the MOI for that cell. If FALSE,
#' the burst size is independent of the MOI.
#' @param choose_strain_by_fitness logical.  If TRUE, strain(s) producing virions
#' in cell is/are chosen proportional to fitness and abundance, otherwise 
#' proportional to abundance only
#' @param one_strain_produced logical. If TRUE, each cell can only produce progeny
#' from one of its co-infecting strains.  If FALSE, each cell can produce progeny
#' for all of its coinfecting strains.
#' @param reassort logical. If TRUE, random reassortment of segments occurs
#'  before packaging.  If FALSE, each strain in the cell effectively produces
#'  pre-packaged viruses
#' @return matrix of dim c(generations, n_strains).  Proportion of virions of
#' each strain for each generation.
#' @importFrom magrittr %>%
#' @export
simulate_evolution <- function(iv, fitness, burst_size, n_cells, 
                               pop_size, generations, mutation_prob,
                               coinfection,
                               MOI_dependent_burst_size,
                               choose_strain_by_fitness,
                               one_strain_produced,
                               reassort) {
  # faster implementation for the case with no co-infection
  if(!coinfection) {
    return(run_mutation_model_no_coinfection(iv, fitness, burst_size, pop_size,
                                             generations, mutation_prob))
  }
  
  n_strains <- 4 # MM, MW, WM, WW
  
  # function to assign cell number to each virion
  assign_cells <- function(virus_popn) {
    cbind(virus_popn, sample.int(n_cells, length(virus_popn), replace=TRUE))
  }
  
  # assign strains to initial virus population and
  # decide which virion is in which cell (assume Poisson distribution)
  virus_popn <- round_preserve_sum(normalise(iv) * pop_size)
  virus_popn <- enumerate_popn(virus_popn)
  virus_popn <- assign_cells(virus_popn)
  
  colnames(virus_popn) <- c("strain", "cell_no")
  
  # initialise results matrix
  results <- matrix(0, nrow = generations, ncol = n_strains)
  results[1,] <- sum_strains(virus_popn[,"strain"])
  
  mutate_popn <- mutate_popn_wrapper(mutation_prob)
  
  # run simulation
  # output_timing <- "timing.txt"
  # old_time <- Sys.time()
  for(generation in 2:generations)  {
    # determine the number of virions produced by a given cell, and what strains they belong to
    make_new_popn <- function(virus_popn) {
      # determine the number of virions of each strain that infected the cell
      strains_in_cell <- sum_strains(virus_popn)
      # determine the burst size from this cell.
      # If the burst size is MOI-independent, then it is burst_size *
      # the average fitness of the co-infecting virions.
      # If the burst size of MOI-dependent, then it is burst_size *
      # the summed fitnesses of the co-infectin virions.
      burst_size_from_cell <- burst_size * sum(fitness * strains_in_cell)
      if(!MOI_dependent_burst_size) {
        burst_size_from_cell <- burst_size_from_cell / sum(strains_in_cell)
      }
      burst_size_from_cell <- rpois(1, burst_size_from_cell)  
      # return early if no virions are produced
      if(burst_size_from_cell == 0) {
        return(numeric(n_strains))
      }
      # for each virion produced, prob_strain[x] gives the probability that 
      # the virion is from strain x
      # if choose_strain_by_fitness is TRUE, prob_strain is a function of both
      # the number of virions of each strain infecting the cell and the fitness
      # of each strain; of FALSE, prob_strain is a function of the number of 
      # virions of each strain infecting the cell only
      if(choose_strain_by_fitness) {
        prob_strain <- normalise(strains_in_cell * fitness)
      } else {
        prob_strain <- normalise(strains_in_cell)
      }
      
      # old implementation of reassortment.  rmultinom at the strain level + reassort
      # produces less variance than rbinom at the segment level
      # if(one_strain_produced) {
      #     new_popn <- as.numeric(rmultinom(1, 1, prob_strain) * burst_size)
      # } else {
      #     new_popn <- as.numeric(rmultinom(1, burst_size, prob_strain))
      # }
      # 
      # # reassort then mutate.  Does the order matter?
      # if(reassort) {
      #     new_popn <- reassort_popn(new_popn)
      # }
      
      # if one_strain_produced is TRUE, choose the strain all the newly produced
      # virions belong to according to prob_strain
      if(one_strain_produced) {
        new_popn <- as.numeric(rmultinom(1, 1, prob_strain) * burst_size_from_cell)
      } else {
        # if one_strain_produced is FALSE and there is reassortment, choose
        # the newly produced segments and randomly package them
        if(reassort) {
          new_popn <- make_and_package_segments(prob_strain, burst_size_from_cell)
        } else {
          # if one_strain_produced is FALSE and there is no reassortment, choose
          # the strain each newly produced virion belongs to according to prob_strain
          new_popn <- as.numeric(rmultinom(1, burst_size_from_cell, prob_strain))
        }
      }
      new_popn
    }
    
    # apply the above function to all infected cells and combine the results
    virus_popn <- tapply(virus_popn[,"strain"], virus_popn[,"cell_no"], make_new_popn)
    virus_popn <- do.call(rbind, virus_popn) %>%
      colSums
    
    # mutate the newly produced strain population
    if(mutation_prob > 0) {
      virus_popn <- mutate_popn(virus_popn)   
    }
    virus_popn <- enumerate_popn(virus_popn)
    
    # cap virus population (assume the effective reproduction number is 1
    # and thus the virus population is constant)
    if(length(virus_popn) > pop_size) {
      virus_popn <- virus_popn[sample.int(length(virus_popn), pop_size)]
    }
    
    # assign which cell each new virion infects.
    # assumes each new generation infects the same number of cells --
    # probably ok if effective reproduction number = 1.
    virus_popn <- assign_cells(virus_popn) 
    colnames(virus_popn) <- c("strain", "cell_no")
    results[generation,] <- sum_strains(virus_popn[,"strain"])
    # )
    # new_time <- Sys.time()
    # timing <- new_time - old_time
    # print(timing)
    # old_time <- new_time
    # write(timing, output_timing, append = TRUE)
  }
  
  # return proportion of each strain at given generation number
  results <- apply(results, 1, normalise)
  rownames(results) <- c("MM", "MW", "WM", "WW")
  t(results)
}

#' make matrix of probabilities of strain x mutating into strain y
#' 
#' @param mutation_prob mutation probability
#' @return 4x4 matrix of mutation probabilities
make_mutation_matrix <- function(mutation_prob) {
  n_strains <- 4
  mutation_matrix <- matrix(c(0, mutation_prob, mutation_prob, mutation_prob^2,
                              mutation_prob, 0, mutation_prob^2, mutation_prob,
                              mutation_prob, mutation_prob^2, 0, mutation_prob,
                              mutation_prob^2, mutation_prob, mutation_prob, 0),
                            nrow = n_strains, byrow = TRUE)
  diag(mutation_matrix) <- 1 - rowSums(mutation_matrix)
  mutation_matrix
}

#' make function to mutate a population of virions
#' 
#' @param mutation_prob probability of mutation for each virion.  Assume same
#' for pb1 and pa.
#' @return function that takes the argument virus_popn: numeric vector of 
#' length 4, with the number of virions in each of the 4 strains, and outputs
#' the same vector but with mutated virions
#' @importFrom magrittr %>%
mutate_popn_wrapper <- function(mutation_prob) {
  mutation_matrix <- make_mutation_matrix(mutation_prob) %>%
    as.data.frame
  function(virus_popn) {
    Map(rmultinom, virus_popn, mutation_matrix, n = 1) %>%
      do.call(cbind, .) %>%
      rowSums
    # round_preserve_sum(as.numeric(mutation_matrix %*% matrix(virus_popn, ncol = 1)))
  }
}

#' old implementation of reassortment
#' 
#' @param virus_popn numeric vector of length 4, with the number of virions
#' in each of the 4 strains
#' @return numeric vector of length 4, with the number of virions
#' in each of the 4 strains, after reassortment

reassort_popn <- function(virus_popn) {
  strain_segments <- matrix(c(1, 1, 1, 0, 0, 1, 0, 0), ncol = 2, byrow = TRUE)
  popn_segments <- enumerate_popn(virus_popn)
  popn_segments <- t(vapply(popn_segments, function(x) strain_segments[x,], numeric(2)))
  popn_segments[,2] <- sample(popn_segments[,2])
  new_popn <- apply(popn_segments, 1, segments_to_strain)
  new_popn <- sum_strains(new_popn)
  new_popn
}

#' take a vector with the number of virions of each strain, and return a vector
#' of length equal to the number of virions, where each entry is the strain number
#' of that virion (the reverse of sum_strains)
#' 
#' @param virus_popn numeric vector of length 4, with the number of virions
#' in each of the 4 strains
#' @return numeric vector of length sum(virus_popn), where each entry is the 
#' strain number of a virion
enumerate_popn <- function(virus_popn) {
  unlist(Map(function(x, y) rep(x, each = y), 
             seq_along(virus_popn),
             virus_popn))
}

#' turns a list of segment indices into a strain number
#' 
#' @param idx numeric vector of length 2, where each entry is 1 or 0: 
#' 1 in the first entry indicates MT PB1,
#' 1 in the second entry indicates MT PA
#' @return the strain number of the virion -- MM = 1, MW = 2, WM = 3, WW = 4
segments_to_strain <- function(idx) {
  -(2 * idx[,1] + idx[,2] + 1) + 5
}

#' takes a vector of length equal to the number of virions, 
#' where each entry is the strain number of that virion, and returns a vector 
#' with the number of virions of each strain (the reverse of enumerate_popn)
#' @param strains a vector of length equal to the number of virions, 
#' where each entry is the strain number of that virion
#' @return numeric vector of length 4, with the number of virions
#' in each of the 4 strains
sum_strains <- function(strains) {
  n_strains <- 4
  vnapply(seq_len(n_strains), function(x) sum(strains == x))
}

#' returns the number of newly produced virions after reassortment and before mutation
#' 
#' @param prob_strain numeric vector of length 4.  Probabilities of the segments
#' corresponding to each of the four strains being produced.
#' @param burst_size burst size from a given cell
#' @return numeric vector of length 4, with the number of newly produced virions
#' in each of the 4 strains
make_and_package_segments <- function(prob_strain, burst_size) {
  strain_segments <- matrix(c(1, 1, 1, 0, 0, 1, 0, 0), ncol = 2, byrow = TRUE)
  prob_wt_segment <- colSums(strain_segments * prob_strain)
  # segments <- vnapply(prob_wt_segment, function(x) rbinom(1, burst_size, x))
  # strains <- numeric(length(prob_strain))
  # # pair WT PA segments[1] up with PB1s -- sample from segments[2] WT PB1s out of burst_size PB1s
  # strains[1] <- rhyper(1, segments[2], burst_size - segments[2], segments[1])
  # # rest of WT PA segments must be paired with mutants
  # strains[3] <- segments[1] - strains[1]
  # remaining_WT_PB1 <- segments[2] - strains[1]
  # remaining_MUT_PB1 <- burst_size - segments[1] - remaining_WT_PB1
  # strains[2] <- rhyper(1, remaining_WT_PB1, remaining_MUT_PB1, burst_size - segments[1])
  # strains[4] <- burst_size - sum(strains)
  # stopifnot(all(strains >= 0))
  
  
  segments <- vapply(prob_wt_segment, function(x) rbinom(burst_size, 1, x), numeric(burst_size))
  if(!is.matrix(segments)) {
    segments <- matrix(segments, nrow = 1)
  }
  strains <- sum_strains(segments_to_strain(segments))
  strains
}