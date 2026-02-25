library(dplyr)
library(ggplot2)

plot_out_dir <- "plots"

if (!dir.exists(plot_out_dir)) {
    dir.create(plot_out_dir)
}

CRS_data <- read.csv("CRS_data.csv")
effect_sizes <- read.csv("all_effect_sizes.csv")

fx_size_species <- effect_sizes %>%
    left_join(CRS_data %>% select(strain, PMID) %>% unique(), by = "PMID", relationship = "many-to-many") %>%
    group_by(Outcome)

pdf(paste0(plot_out_dir, "/effect_size_by_species.pdf"), width = 12, height = 8)

ggplot(
    fx_size_species,
    aes(x = strain, y = effect_size, fill = strain)
) +
    geom_boxplot() +
    facet_wrap(~Outcome, scales = "free_x") +
    theme_bw() +
    ylim(-10, 10) +
    xlab("") +
    ylab("Effect size (Hedges' g)") +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 15, face = "bold")
    )

dev.off()

# p = 0.80
aov(effect_size ~ strain, data = fx_size_species %>% 
    filter(Outcome == "EPM")) %>% 
    summary()

# p = 0.83
lm(effect_size ~ strain, data = fx_size_species %>% 
    filter(Outcome == "FST")) %>% 
    summary()

# p = 0.23
lm(effect_size ~ strain, data = fx_size_species %>% 
    filter(Outcome == "OFT")) %>% 
    summary()

# = 0.23
lm(effect_size ~ strain, data = fx_size_species %>% 
    filter(Outcome == "SPT")) %>% 
    summary()
