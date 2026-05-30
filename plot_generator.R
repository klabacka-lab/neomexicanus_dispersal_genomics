library(tidyverse)
library(tidyr)
library(patchwork)
library(ggplot2)

data <- read.csv("/home/perry/Desktop/Classes/klab/local_code/neomex/avg_kmer_filesize.csv", header = FALSE)
colnames(data) <- c("kmer", "file_size")

p1 <- ggplot(data, aes(x = factor(kmer), y = file_size)) +
  geom_col() +
  coord_cartesian(ylim = c(10.6, 11.0)) +
  labs(
    x = NULL,
    y = "Average File Size (GB)",
    title = "Xengsort Output size by K-mer Length"
  ) +
  theme_minimal()

data$delta_size <- c(NA, diff(data$file_size))

p2 <- ggplot(data[-1, ], aes(x = kmer, y = delta_size)) +
  geom_col() +
  labs(
    x = "k-mer size",
    y= "File Size Increase",
    title = "Average File Size increase From Previous K"
  ) +
  theme_minimal()

p1 / p2
