library(latex2exp)
library(tidyverse)

dir_name <- "~/git_repos/reassortment/results_0point01/sensitivity/"
dir_name_plot <- "~/overleaf/reassortment/sensitivity/figs_revision1/"
default_w_reassortment_filename <- paste0(dir_name, "default/results_A.rds")
default_wo_reassortment_filename <- paste0(dir_name, "D/results_D.rds")
fitness_WM <- c(1, 1.1, 1.25, 1.5, 1.75, 2)
reassort <- c(TRUE, FALSE)
par_grid <- expand.grid(fitness_WM = fitness_WM, reassort = reassort)
filenames <- paste0(dir_name, "fitness_WM_", num2str(par_grid$fitness_WM), "_",
                    par_grid$reassort, "/results.rds")
# didn't re-run defaults, so the default filenames have a different naming system
filenames[par_grid$fitness_WM == 1.25 & par_grid$reassort] <- default_w_reassortment_filename
filenames[par_grid$fitness_WM == 1.25 & !(par_grid$reassort)] <- default_wo_reassortment_filename

results <- lapply(filenames, readRDS)
g <- sensitivity_plot(results, par_grid$fitness_WM, par_grid$reassort, x_lab = TeX("$f_{WM}$"))
filename <- paste0(dir_name_plot, "f_WM.pdf")
ggsave(filename, g$g, width = 15, height = 10, units = "cm")

fitness_MM <- c(0.8, 0.9, 1, 1.1, 1.2)
par_grid <- expand.grid(fitness_MM = fitness_MM, reassort = reassort)
filenames <- paste0(dir_name, "fitness_MM_", num2str(par_grid$fitness_MM), "_",
                    par_grid$reassort, "/results.rds")
filenames[par_grid$fitness_MM == 1 & par_grid$reassort] <- default_w_reassortment_filename
filenames[par_grid$fitness_MM == 1 & !(par_grid$reassort)] <- default_wo_reassortment_filename
results <- lapply(filenames, readRDS)
g <- sensitivity_plot(results, par_grid$fitness_MM, par_grid$reassort, x_lab = TeX("$f_{MM}$"))
filename <- paste0(dir_name_plot, "f_MM.pdf")
ggsave(filename, g$g, width = 15, height = 10, units = "cm")

mu <- c(2e-3, 2e-4, 2e-5)
par_grid <- expand.grid(mu = mu, reassort = reassort)
filenames <- paste0(dir_name, "mutation_prob_", num2str(par_grid$mu), "_",
                    par_grid$reassort, "/results.rds")
filenames[par_grid$mu == 2e-4 & par_grid$reassort] <- default_w_reassortment_filename
filenames[par_grid$mu == 2e-4 & !(par_grid$reassort)] <- default_wo_reassortment_filename
results <- lapply(filenames, readRDS)
g <- sensitivity_plot(results, par_grid$mu, par_grid$reassort, x_lab = "Mutation rate per generation")
filename <- paste0(dir_name_plot, "mu.pdf")
ggsave(filename, g$g, width = 15, height = 10, units = "cm")

sim_name <- c("default","MOI_dependent_burst_size_FALSE", "choose_strain_by_fitness_TRUE")
sim_name_plot <- c("default", "burst size\n independent of MOI", "progeny composition\n depends on fitness")
par_grid <- expand.grid(sim_name = sim_name, reassort = reassort)
filenames <- paste0(dir_name, sim_name, "_",
                    par_grid$reassort, "/results.rds")
filenames[par_grid$sim_name == "default" & par_grid$reassort] <- default_w_reassortment_filename
filenames[par_grid$sim_name == "default" & !(par_grid$reassort)] <- default_wo_reassortment_filename
results <- lapply(filenames, readRDS)
par_grid <- expand.grid(sim_name = sim_name_plot, reassort = reassort)
g <- sensitivity_plot(results, par_grid$sim_name, par_grid$reassort, x_lab = "", rotate_x_labels = TRUE)
filename <- paste0(dir_name_plot, "misc.pdf")
ggsave(filename, g$g, width = 15, height = 15, units = "cm")

fitness_MW <- c(.01, .1, .5, 1)
par_grid <- expand.grid(fitness_MW = fitness_MW, reassort = reassort)
filenames <- paste0(dir_name, "fitness_MW", (par_grid$fitness_MW), 
                    par_grid$reassort, "/results.rds")
filenames[par_grid$fitness_MW == .01 & par_grid$reassort] <- default_w_reassortment_filename
filenames[par_grid$fitness_MW == .01 & !(par_grid$reassort)] <- default_wo_reassortment_filename
results <- lapply(filenames, readRDS)
g <- sensitivity_plot(results, par_grid$fitness_MW, par_grid$reassort, x_lab = "Fitness of PB1 K229R")
filename <- paste0(dir_name_plot, "fitness_MW.pdf")
ggsave(filename, g$g, width = 15, height = 10, units = "cm")

iv_WT_double_mutant <- c(seq(5, 50, by = 10))
par_grid <- expand.grid(iv_WT_double_mutant = iv_WT_double_mutant, reassort = reassort)
filenames <- paste0(dir_name, "iv_WT_double_mutant_", (par_grid$iv_WT_double_mutant), "_",
                    100 - par_grid$iv_WT_double_mutant,
                    par_grid$reassort, "/results.rds")
filenames[par_grid$iv_WT_double_mutant == 5 & par_grid$reassort] <- default_w_reassortment_filename
filenames[par_grid$iv_WT_double_mutant == 5 & !(par_grid$reassort)] <- default_wo_reassortment_filename
results <- lapply(filenames, readRDS)
g <- sensitivity_plot(results, par_grid$iv_WT_double_mutant/100, par_grid$reassort, x_lab = "Initial proportion of WT")
filename <- paste0(dir_name_plot, "iv_WT_double_mutant.pdf")
ggsave(filename, g$g, width = 15, height = 10, units = "cm")

prop_WT <- seq(5, 15, by = 5)
par_grid <- expand.grid(prop_WT = prop_WT, reassort = reassort)
filenames <- paste0(dir_name, "iv_", 100 - par_grid$prop_WT, "_", 
                    par_grid$prop_WT, "_", 
                    par_grid$prop_WT, "_", 
                    par_grid$prop_WT, "_",
                    par_grid$reassort, "/results.rds")
results <- lapply(filenames, readRDS)
g <- sensitivity_plot(results, par_grid$prop_WT/100, par_grid$reassort, x_lab = "Initial proportion of WT, PB1 K229R, PA P653L")
filename <- paste0(dir_name_plot, "prop_WT.pdf")
ggsave(filename, g$g, width = 15, height = 10, units = "cm")

MOI <- c(.1, 1, 10)
par_grid <- expand.grid(MOI = MOI, reassort = reassort)
filenames <- paste0(dir_name, "MOI_", num2str(par_grid$MOI), "_",
                    par_grid$reassort, "/results.rds")
filenames[par_grid$MOI == 1 & par_grid$reassort] <- default_w_reassortment_filename
filenames[par_grid$MOI == 1 & !(par_grid$reassort)] <- default_wo_reassortment_filename
results <- lapply(filenames, readRDS)
g <- sensitivity_plot(results, par_grid$MOI, par_grid$reassort, x_lab = "MOI")
filename <- paste0(dir_name_plot, "MOI.pdf")
ggsave(filename, g$g, width = 15, height = 10, units = "cm")

pop_size <- c(1e5, 1e6, 5e6)
par_grid <- expand.grid(pop_size = pop_size, reassort = reassort)
filenames <- paste0(dir_name, "pop_size", (par_grid$pop_size), par_grid$reassort, "/results.rds")
filenames[par_grid$pop_size == 1e6 & par_grid$reassort] <- default_w_reassortment_filename
filenames[par_grid$pop_size == 1e6 & !(par_grid$reassort)] <- default_wo_reassortment_filename
results <- lapply(filenames, readRDS)
g <- sensitivity_plot(results, par_grid$pop_size, par_grid$reassort, x_lab = "Viral load (virions)")
filename <- paste0(dir_name_plot, "pop_size.pdf")
ggsave(filename, g$g, width = 15, height = 10, units = "cm")
