library(dplyr)
library(robvis)

SYRCLE <- read.csv("SYRCLE analysis.csv")

# We need to format as expected in robvis

SYRCLE %>% 
    mutate(FirstAuthor = unlist(lapply(sapply(SYRCLE$Authors, strsplit, " "), "[", 1))) %>%
    mutate(Study = paste(FirstAuthor, Publication.Year, sep = " ")) %>%
    # Rename questions to D1, D2, ..., D10
    rename_with(~ paste0("D", 1:10), 7:16) %>%
    select(Study, starts_with("D"), -DOI) %>%
    # Change Yes, Data complete or Similar at baseline to Low, No to High
    mutate(across(starts_with("D"), ~ case_when(
        . %in% c("Yes", "Data complete", "Similar at baseline") ~ "Low",
        . == "No" ~ "High",
        TRUE ~ "Unclear"
    ))) %>% 
    rob_traffic_light(tool = "ROB2")
