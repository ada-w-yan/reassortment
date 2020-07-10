library(latex2exp)
dir_name <- "~/overleaf/reassortment/sensitivity/figs/"
fitness_WM <- c(1, 1.1, 1.5, 1.75, 2)
reassort <- c(TRUE, FALSE)
par_grid <- expand.grid(fitness_WM = fitness_WM, reassort = reassort)
filenames <- paste0(dir_name, "fitness_WM_", num2str(par_grid$fitness_WM), "_",
                    par_grid$reassort, "/results.rds")
results <- lapply(filenames, readRDS)
g <- sensitivity_plot(results, par_grid$fitness_WM, par_grid$reassort, x_lab = TeX("$f_{WM}$"))
filename <- paste0(dir_name, "f_WM.pdf")
ggsave(filename, g$g, width = 15, height = 10, units = "cm")

fitness_MM <- c(0.8, 0.9, 1.1, 1.2)
par_grid <- expand.grid(fitness_MM = fitness_MM, reassort = reassort)
filenames <- paste0(dir_name, "fitness_MM_", num2str(par_grid$fitness_MM), "_",
                    par_grid$reassort, "/results.rds")
results <- lapply(filenames, readRDS)
g <- sensitivity_plot(results, par_grid$fitness_MM, par_grid$reassort, x_lab = TeX("$f_{MM}$"))
filename <- paste0(dir_name, "f_MM.pdf")
ggsave(filename, g$g, width = 15, height = 10, units = "cm")

mu <- c(2e-3, 2e-5)
par_grid <- expand.grid(mu = mu, reassort = reassort)
filenames <- paste0(dir_name, "mutation_prob_", num2str(par_grid$mu), "_",
                    par_grid$reassort, "/results.rds")
results <- lapply(filenames, readRDS)
g <- sensitivity_plot(results, par_grid$mu, par_grid$reassort, x_lab = TeX("Mutation rate (day$^{-1}$)"))
filename <- paste0(dir_name, "mu.pdf")
ggsave(filename, g$g, width = 15, height = 10, units = "cm")

sim_name <- c("MOI_0point1", "MOI_dependent_burst_size_FALSE", "choose_strain_by_fitness_TRUE")
sim_name_plot <- c("MOI = 0.1", "MOI independent of\n burst size", "progeny composition\n depends on fitness")
par_grid <- expand.grid(sim_name = sim_name, reassort = reassort)
filenames <- paste0(dir_name, sim_name, "_",
                    par_grid$reassort, "/results.rds")
results <- lapply(filenames, readRDS)
par_grid <- expand.grid(sim_name = sim_name_plot, reassort = reassort)
g <- sensitivity_plot(results, par_grid$sim_name, par_grid$reassort, x_lab = "", rotate_x_labels = TRUE)
filename <- paste0(dir_name, "misc.pdf")
ggsave(filename, g$g, width = 15, height = 15, units = "cm")
