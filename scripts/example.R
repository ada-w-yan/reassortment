library(magrittr)
library(ggplot2)
set.seed(2)
iv <- c(95,0,0,5) #mt,mt  mt,wt  wt,mt   wt,wt
fitness <- c(1,0,1.1,1)
n_cells <- 1e4
pop_size <- 1000
mutation_prob <- 0#1e-3
burst_size <- 10
generations <- 40

# results <- run_reassortment_model(iv, fitness, burst_size, n_cells, pop_size, generations)
# plot_strains(results)
# plot_segments(results)

# # results <- run_mutation_model_no_coinfection(iv, fitness, burst_size, pop_size,
#                                               generations, mutation_prob)
results <- run_mutation_model(iv, fitness, burst_size, n_cells, pop_size,
                                             generations, mutation_prob,
                                            choose_strain_by_fitness = FALSE,
                              one_strain_produced = FALSE,
                              reassort = TRUE)

# results are probably dependent on:
# MOI = population size / n_cells
# generations

# results are probably not dependent on:
# burst_size
# n_cells
# fitness
# mutation model: disallow co-infection? then how to calibrate the model to have the same average burst size?
# set burst size in new model to be burst size in old model * mean MOI?
# Subsequent generation only depends on one virion inside the infected cell (which boils down to the same thing)?
# subseqent generation burst size depends on all strains, but only one strain can be produced from the infected cell?
# subsequent generation burst size depends on all strains, and a co-infected cell can produce both strains but not reassorted strains?
