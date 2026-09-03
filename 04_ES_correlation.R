library(dplyr)
library(tidyr)
library(corrplot)

out_dir <- "plots"

CRS_effects <- read.csv("all_effect_sizes.csv")

pdf(file.path(out_dir, "correlation_plot.pdf"), 
    width = 10, height = 10)

# Ampuero 2015 would result as a duplicate, because the difference is in how they restrain the animals
# We give them different IDs
CRS_effects[CRS_effects$ID == 9663, "ID"] <- c(96631, 96632, 96631, 96632)
CRS_effects %>% 
    select(ID, Outcome, effect_size, total_time_h, sex, strain) %>%
    pivot_wider(
        id_cols = c(ID, total_time_h, sex, strain),
        names_from = Outcome, values_from = effect_size
    ) %>%
    select(-c(ID, total_time_h, sex, strain)) %>%
    cor(use = "pairwise.complete.obs", method = "pearson") %>%
    corrplot.mixed(
        upper = "ellipse", lower = "number",
        order = "hclust", upper.col = COL2("PRGn"),
        lower.col = "black", tl.col = "black",
        tl.cex = 2, cl.cex = 1.2, number.cex = 1.6
    )

dev.off()
