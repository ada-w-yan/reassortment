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
#' @importFrom dplyr %>%
#' @export
plot_multirun_strains <- function(results) {
  colours <- c("purple", "red", "blue", "black")
  results <- as.data.frame(results)
  names(colours) <- colnames(results)[-ncol(results)]
  n_runs <- results[nrow(results), "run"]
  n_gen <- nrow(results) / n_runs
  results$gen <- seq_len(n_gen)
  g <- results %>% 
    tidyr::gather(key = "strain", value = "proportion", -gen, -run) %>%
    dplyr::group_by(gen, strain) %>%
    dplyr::summarise(quantile(proportion, 0.025), quantile(proportion, 0.975)) %>%
    ggplot(aes(x = gen, 
               ymin = `quantile(proportion, 0.025)`,
               ymax = `quantile(proportion, 0.975)`,
               color = strain, 
               group = strain)) +
    geom_errorbar() + 
    theme_bw() +
    coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    xlab("Generation") +
    ylab("Proportion of strain") +
    scale_color_manual("", values = colours) +
    theme(legend.justification=c(1,1),
          legend.position=c(.9,.5),
          text = element_text(size = 20))
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