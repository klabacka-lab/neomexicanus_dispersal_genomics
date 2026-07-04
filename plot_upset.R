library(ComplexUpset)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

input <- args[1]
output <- args[2]

df <- read.csv(input, check.names = FALSE)
rownames(df) <- NULL

sample_cols <- setdiff(colnames(df), c("CHROM","POS"))

pdf(output, width=10, height=6)

upset(
      df,
      intersect = sample_cols,
      name = "Individuals",
      width_ratio = 0.25
)

dev.off()
