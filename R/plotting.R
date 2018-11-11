#' plot results of run_xxx_model by strain
#' 
#' @param results output of run_xxx_model  
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

#' plot results of run_xxx_model by segment
#' 
#' @param results output of run_xxx_model  
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