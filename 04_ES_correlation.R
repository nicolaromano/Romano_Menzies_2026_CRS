library(dplyr)
library(tidyr)
library(corrplot)

out_dir <- "plots"

CRS_effects <- read.csv("all_effect_sizes.csv")

pdf(file.path(out_dir, "correlation_plot.pdf"), 
    width = 10, height = 10)

CRS_effects %>%
    select(PMID, Outcome, effect_size, total_time_h, sex, strain) %>%
    pivot_wider(
        id_cols = c(PMID, total_time_h, sex, strain),
        names_from = Outcome, values_from = effect_size
    ) %>%
    select(-c(PMID, total_time_h, sex, strain)) %>%
    cor(use = "pairwise.complete.obs", method = "pearson") %>%
    corrplot.mixed(
        upper = "ellipse", lower = "number",
        order = "hclust", upper.col = COL2("PRGn"),
        lower.col = "black", tl.col = "black",
        tl.cex = 2, cl.cex = 1.2, number.cex = 1.6
    )

dev.off()
