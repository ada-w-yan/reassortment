#' plot results of run_xxx_model by strain
#' 
#' @param results output of run_xxx_model  
#' @import ggplot2
#' @export
plot_strains <- function(results, zoom = FALSE) {
  colours <- c("purple", "red", "blue", "black")
  results <- as.data.frame(results)
  results$gen <- seq_len(nrow(results))
  results_melt <- reshape2::melt(results, id.var = "gen")
  ymax <- ifelse(zoom, max(results$WM), 1)
  ggplot(results_melt, aes(x = gen, y = value, color = variable, group = variable)) +
    geom_line() +
    theme_bw() +
    coord_cartesian(ylim = c(0, ymax), expand = FALSE) +
    xlab("Generation") +
    ylab("Proportion of viral load") +
    scale_colour_manual("Strain", values = colours) +
    theme(legend.justification=c(1,1),
          legend.position=c(.9,.5),
          text = element_text(size = 20))
}

#' plot results of run_xxx_model by segment
#' 
#' @param results output of run_xxx_model  
#' @import ggplot2
#' @export
plot_segments <- function(results) {
  colours <- c("red", "blue")
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
    ylab("Proportion of viral load") +
    scale_colour_manual("Segment", labels = c("MT PB1", "MT PA"), values = colours) +
    theme(legend.justification=c(1,1),
          legend.position=c(.9,.5),
          text = element_text(size = 20))
  
}