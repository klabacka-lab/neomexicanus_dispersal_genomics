library(ComplexUpset)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

input <- args[1]
output <- args[2]

df <- read.csv(input, check.names = FALSE)
rownames(df) <- NULL

sample_cols <- colnames(df)[-1]
sample_cols <- c("7888", "7890", "7891", "7892", "7893", "7894", "7895", "7896", "7897")

pdf(output, width=9, height=7)

upset(
      df,
      intersect = sample_cols,
      name = "Individuals",
      width_ratio = 0.25,
      n_intersections = 15,
      min_degree = 2,
      sort_intersections_by = "cardinality"
)

dev.off()
