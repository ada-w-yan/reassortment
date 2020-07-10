#' plot results of run_xxx_model by strain
#' 
#' @param results output of run_xxx_model  
#' @import ggplot2
#' @importFrom dplyr %>%
#' @export
plot_strains <- function(results, zoom = FALSE) {
  colours <- c("purple", "red", "blue", "black")
  results <- as.data.frame(results)
  names(colours) <- colnames(results)
  results$gen <- seq_len(nrow(results))
  ymax <- ifelse(zoom, max(results$WM), 1)
  g <- results %>% 
    tidyr::gather(key = "strain", value = "proportion", -gen) %>%
  ggplot(aes(x = gen, y = proportion, color = strain, group = strain)) +
    geom_line() +
    theme_bw() +
    coord_cartesian(ylim = c(0, ymax), expand = FALSE) +
    xlab("Generation") +
    ylab("Proportion of viral load") +
    scale_colour_manual("", values = colours) +
    theme(legend.justification=c(1,1),
          legend.position=c(.9,.5),
          text = element_text(size = 20))
}

#' plot results of run_xxx_model by strain when there are multiple stochastic runs
#' 
#' @param results output of run_xxx_model  
#' @import ggplot2
#' @importFrom dplyr group_by summarise
#' @importFrom tidyr pivot_longer
#' @export
plot_multirun_strains <- function(results) {
  colours <- c("purple", "red", "blue", "black")
  results <- as.data.frame(results)
  names(colours) <- colnames(results)[-ncol(results)]
  n_runs <- results[nrow(results), "run"]
  n_gen <- nrow(results) / n_runs
  results$gen <- seq_len(n_gen)
  g <- results %>% 
    pivot_longer(c(-gen, -run), names_to = "strain", values_to = "proportion") %>%
    group_by(gen, strain) %>%
    summarise(quantile(proportion, 0.025), quantile(proportion, 0.975)) %>%
    ggplot(aes(x = gen, 
               ymin = `quantile(proportion, 0.025)`,
               ymax = `quantile(proportion, 0.975)`,
               color = strain,
               fill = strain, 
               group = strain)) +
    geom_ribbon() + 
    theme_bw() +
    scale_x_continuous(breaks = c(1, 5, 10, 15, 20), limits = c(0, 20), expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1), expand = c(0,0)) +
    xlab("Time (generations)") +
    ylab("Proportion") +
    scale_fill_manual("", values = colours) +
    scale_color_manual("", values = colours) +
    theme(legend.position = "none",
    # theme(legend.justification=c(1,1),
    #       legend.position=c(.9,.5),
          text = element_text(size = 16),
    plot.margin = margin(.5, .5, .5, .5, "cm"))
    
  g
}

#' plot results of run_xxx_model by segment
#' 
#' @param results output of run_xxx_model  
#' @import ggplot2
#' @importFrom dplyr %>%
#' @export
plot_segments <- function(results) {
  colours <- c(PB1 = "red", PA = "blue")
  results <- as.data.frame(results)
  results$gen <- seq_len(nrow(results))
  g <- results %>% 
    dplyr::mutate(PB1 = MM + MW, PA = MM + WM) %>%
    tidyr::gather(key = "segment", value = "proportion", PB1, PA) %>%
  ggplot(aes(x = gen, y = proportion, color = segment, group = segment)) +
    geom_line() +
    theme_bw() +
    coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    xlab("Generation") +
    ylab("Proportion of mutant segment") +
    scale_colour_manual("", values = colours) +
    theme(legend.justification=c(1,1),
          legend.position=c(.9,.5),
          text = element_text(size = 20))
  
}

#' plot results of run_xxx_model by segment when there are multiple stochastic runs
#' 
#' @param results output of run_xxx_model  
#' @import ggplot2
#' @importFrom dplyr %>%
#' @export
plot_multirun_segments <- function(results) {
  colours <- c(PB1 = "red", PA = "blue")
  results <- as.data.frame(results)
  n_runs <- results[nrow(results), "run"]
  n_gen <- nrow(results) / n_runs
  results$gen <- seq_len(n_gen)
  g <- results %>% 
    dplyr::mutate(PB1 = MM + MW, PA = MM + WM) %>%
    tidyr::gather(key = "segment", value = "proportion", PB1, PA) %>%
    dplyr::group_by(gen, segment) %>%
    dplyr::summarise(quantile(proportion, 0.025), quantile(proportion, 0.975)) %>%
  ggplot(aes(x = gen, 
             ymin = `quantile(proportion, 0.025)`,
             ymax = `quantile(proportion, 0.975)`,
             color = segment, 
             group = segment)) +
    geom_errorbar() + 
    theme_bw() +
    coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    xlab("Generation") +
    ylab("Proportion of mutant segment") +
    scale_color_manual("", values = colours) +
    theme(legend.justification=c(1,1),
          legend.position=c(.9,.5),
          text = element_text(size = 20))
}

#' plot proportion of each strain at 20 generations across sensitivity analyses
#' 
#' @param results_list list of results; each element is an output of run_xxx_model
#' @param results_name vector of names of simulation to plot along x-axis
#' @param reassortment_in vector of logicals: whether those results were from a simulation
#' with or without reassortment
#' @import ggplot2
#' @importFrom dplyr %>%
#' @export
sensitivity_plot <- function(results_list, results_name, reassortment_in, x_lab, rotate_x_labels = FALSE) {
  colours <- c("purple", "red", "blue", "black")
  collate_results <- function(results, results_name, reassortment_in) {
    n_runs <- results[nrow(results), "run"]
    n_gen <- nrow(results) / n_runs
    results <- as.data.frame(results)
    results$gen <- seq_len(n_gen)
    results %>% as_tibble %>%
      filter(gen == n_gen) %>%
      summarise_all(median) %>%
      select(-run, -gen) %>%
      mutate(reassortment = reassortment_in, sim_name = results_name)
  }
  results <- Map(collate_results, results_list, results_name, reassortment_in) %>%
    do.call(rbind, .)
  names(colours) <- colnames(results)[seq_along(colours)]
  results <- results %>%
    pivot_longer(-c("sim_name", "reassortment")) %>%
    mutate(sim_name = factor(sim_name, levels = results_name[seq_along(unique(results_name))]),
           reassortment = factor(reassortment))
  levels(results$reassortment) <- c("Mutation only", "Reassortment and mutation")
  g <- ggplot(results, aes(x = sim_name, y = value, fill = name)) +
    facet_wrap(~reassortment, nrow = 1) +
    geom_bar(position = "stack", stat = "identity") +
    theme_bw() +
    xlab(x_lab) +
    ylab("Proportion") +
    scale_fill_manual("", values = colours) +
    scale_color_manual("", values = colours) +
    theme(legend.position = "none",
          # theme(legend.justification=c(1,1),
          #       legend.position=c(.9,.5),
          text = element_text(size = 16),
          plot.margin = margin(.5, .5, .5, .5, "cm"))
  if(rotate_x_labels) {
    g <- g + theme(axis.text.x = element_text(angle = 90))
  }
  
  list(g = g, results = results)
}