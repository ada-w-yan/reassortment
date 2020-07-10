wrapper <- function(dir_name) {
  results <- readRDS(paste0(dir_name, "/results.rds"))
  g <- plot_multirun_strains(results)
  ggsave(paste0(dir_name, "/strains.pdf"), g, width = 10, height = 10, units = "cm")
  saveRDS(g, paste0(dir_name, "/strains.rds"))
  g <- plot_multirun_segments(results)
  ggsave(paste0(dir_name, "/segments.pdf"), g, width = 10, height = 10, units = "cm")
  saveRDS(g, paste0(dir_name, "/segments.rds"))
  invisible(g)
}
