#' Simulate evolution of mixture by reassortment
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
#' @return matrix of dim c(generations, n_strains).  Proportion of virions of
#' each strain for each generation.
#' @importFrom magrittr %>%
#' @export
run_reassortment_model <- function(iv, fitness, burst_size, n_cells, pop_size, generations) {
    n_strains <- 4
    # assign strains to initial virus population
    strain <- sample.int(n_strains, 
                         size=pop_size, 
                         replace=TRUE, prob=normalise(iv))
    strain_segments <- matrix(c(1, 1, 1, 0, 0, 1, 0, 0), ncol = 2, byrow = TRUE)
    assign_cells <- function(virus_popn) {
        cbind(virus_popn, sample.int(n_cells, nrow(virus_popn), replace=TRUE))
    }
    # decide which virion is in which cell (assume Poisson distribution)
    virus_popn <- strain_segments[strain,] %>%
        assign_cells
    colnames(virus_popn) <- c("mt_pb1", "mt_pa", "cell_no")
    
    
    
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
            do.call(rbind, .)
        
        # cap virus population at maximum population
        # if effective reproduction number = 1, virus population is constant -- set
        # initial population equal to max population
        if(nrow(virus_popn) > pop_size) {
            virus_popn <- virus_popn[sample.int(nrow(virus_popn), pop_size),]
        }
        
        virus_popn <- assign_cells(virus_popn) # assumes each new generation infects the same number of cells --
        # probably ok if effective reproduction number = 1
        colnames(virus_popn) <- c("mt_pb1", "mt_pa", "cell_no")
        
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
#' @importFrom magrittr %>%
#' @export
run_mutation_model_no_coinfection <- function(iv, fitness, burst_size, pop_size,
                                              generations, mutation_prob) {
    # equivalent to allowing co-infection, but only one virion is selected to be
    # productive and that virion dictates the burst size
    n_strains <- 4
    # assign strains to initial virus population
    popn_by_strain <- as.numeric(rmultinom(1, size=pop_size, prob=normalise(iv)))
    burst_size_by_strain_per_virus <- burst_size * fitness
    results <- matrix(0, nrow = generations, ncol = n_strains)
    results[1,] <- popn_by_strain
    
    mutation_matrix <- make_mutation_matrix(mutation_prob)
    
    for(generation in 2:generations)  {
        burst_size <- probabilistic_round(burst_size_by_strain_per_virus * popn_by_strain) %>%
            matrix(ncol = 1)
        # deterministic mutation
        burst_size <- (mutation_matrix %*% burst_size) %>%
            as.numeric
        # should be selecting from, rather than sampling with probability...
        popn_by_strain <- rmultinom(1, pop_size, normalise(burst_size)) %>%
            as.numeric
        results[generation,] <- popn_by_strain
    }
    
    results <- apply(results, 1, normalise)
    
    rownames(results) <- c("MM", "MW", "WM", "WW")
    t(results)
}

#' Simulate evolution of mixture by mutation, where all progeny by a cell is of
#' the same strain, and that strain is selected randomly from the intracellular
#' population
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
#' @param choose_strain_by_fitness logical.  If TRUE, strain(s) producing virions
#' in cell is/are chosen proportional to fitness and abundance, otherwise 
#' proportional to abundance only
#' @param one_strain_produced logical. If TRUE,
#' @return matrix of dim c(generations, n_strains).  Proportion of virions of
#' each strain for each generation.
#' @importFrom magrittr %>%
#' @export
run_mutation_model <- function(iv, fitness, burst_size, n_cells, 
                               pop_size, generations, mutation_prob,
                               choose_strain_by_fitness,
                               one_strain_produced,
                               reassort) {
    n_strains <- 4
    
    assign_cells <- function(virus_popn) {
        cbind(virus_popn, sample.int(n_cells, length(virus_popn), replace=TRUE))
    }
    
    # assign strains to initial virus population and
    # decide which virion is in which cell (assume Poisson distribution)
    
    virus_popn <- sample.int(n_strains, 
                             size=pop_size, 
                             replace=TRUE, prob=normalise(iv)) %>%
        assign_cells
    
    colnames(virus_popn) <- c("strain", "cell_no")
    
    results <- matrix(0, nrow = generations, ncol = n_strains)
    results[1,] <- sum_strains(virus_popn[,"strain"])
    for(generation in 2:generations)  {
        infected_cells <- unique(virus_popn[,"cell_no"])
        make_new_popn <- function(infected_cell) {
            virions_in_cell <- virus_popn[,"cell_no"]==infected_cell
            strains_in_cell <- virus_popn[virions_in_cell,"strain"] %>%
                sum_strains
            # assume actual burst size is a function of average burst size and fitnesses of each virus in the cell
            # and is deterministic
            burst_size_from_cell <- probabilistic_round(burst_size*sum(fitness * strains_in_cell))  
            if(burst_size_from_cell == 0) {
                return(numeric(n_strains))
            }
            if(choose_strain_by_fitness) {
                prob_strain <- normalise(strains_in_cell * fitness)
            } else {
                prob_strain <- normalise(strains_in_cell)
            }
            if(one_strain_produced) {
                new_popn <- as.numeric(rmultinom(1, 1, prob_strain) * burst_size)
            } else {
                new_popn <- as.numeric(rmultinom(1, burst_size, prob_strain))
            }
            
            # reassort then mutate.  Does the order matter?
            if(reassort) {
                new_popn <- reassort_popn(new_popn)
            }
            new_popn <- mutate_popn(new_popn, mutation_prob)
            new_popn
        }
        
        virus_popn <- vapply(infected_cells, make_new_popn, numeric(n_strains)) %>%
            rowSums %>%
            enumerate_popn
        # cap virus population at maximum population
        # if effective reproduction number = 1, virus population is constant -- set
        # initial population equal to max population
        if(length(virus_popn) > pop_size) {
            virus_popn <- virus_popn[sample.int(length(virus_popn), pop_size)]
        }
        
        virus_popn <- assign_cells(virus_popn) # assumes each new generation infects the same number of cells --
        # probably ok if effective reproduction number = 1
        colnames(virus_popn) <- c("strain", "cell_no")
        results[generation,] <- sum_strains(virus_popn[,"strain"])
    }
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

mutate_popn <- function(virus_popn, mutation_prob) {
    mutation_matrix <- make_mutation_matrix(mutation_prob)
    (mutation_matrix %*% matrix(virus_popn, ncol = 1)) %>%
        as.numeric %>%
        round_preserve_sum
}

reassort_popn <- function(virus_popn) {
    strain_segments <- matrix(c(1, 1, 1, 0, 0, 1, 0, 0), ncol = 2, byrow = TRUE)
    popn_segments <- enumerate_popn(virus_popn) %>%
        vapply(., function(x) strain_segments[x,], numeric(2)) %>%
        t
    popn_segments[,2] <- sample(popn_segments[,2])
    new_popn <- apply(popn_segments, 1, segments_to_strain) %>%
        sum_strains
    new_popn
}

enumerate_popn <- function(virus_popn) {
    unlist(Map(function(x, y) rep(x, each = y), 
               seq_along(virus_popn),
               virus_popn))
}

segments_to_strain <- function(idx) {
    -(2 * idx[1] + idx[2] + 1) + 5
}

sum_strains <- function(strains) {
    vnapply(seq_len(n_strains), function(x) sum(strains == x))
}