library(xlsx)
library(dplyr)

get_data <- function(test, extra_cols=NULL) {    
    data_original <- read.xlsx("08 csv-chronicres-set - test characteristics and effect sizes.xlsx", 
        sheetName = test, startRow = 2)    

    res <- data_original %>% 
        # Rename columns for consistency
        rename(Journal = "Journal.Book",
            Date = "Publication.Year",
            age = "initial.age..weeks.",
            control_mean = "control.mean",
            control_sd = "back.calculated.SD",
            control_n = "n.used",
            test_mean = "test.mean",
            test_sd = "back.calculated.SD.1",
            test_n = "n.used.1",
            total_time = "total.number.of.min",
            total_days = "total.number.of.days"
            ) %>%
        # rename extra_cols to extra_cols_names
        rename(!!!extra_cols) %>%
        # Name publications as "Author Year", so it's not too long
        mutate(First_Author = sub(",.*", "", Authors)) %>%
        mutate(First_Author = sub(" .*", "", First_Author)) %>%
        mutate(PubID = paste(First_Author, Date)) %>%
        mutate(sex = sub("not given", NA, sex)) %>% 
        mutate(control_n = ifelse(is.na(control_n), n, control_n)) %>%
        mutate(test_n = ifelse(is.na(test_n), n.1, test_n)) %>%        
        mutate(control_mean = ifelse(grepl("[a-zA-Z]", control_mean), NA, control_mean)) %>%        
        mutate(control_mean = as.numeric(control_mean)) %>%
        mutate(test_mean = ifelse(grepl("[a-zA-Z]", test_mean), NA, test_mean)) %>%
        mutate(test_mean = as.numeric(test_mean)) %>%
        mutate(Outcome = test) %>%
        mutate(total_days = as.numeric(total_days),
               total_time = as.numeric(total_time)) %>%
        select(PMID, PubID, Journal, Outcome,
                strain, sex, age, 
                control_mean, control_sd, control_n, 
                test_mean, test_sd, test_n, total_days, total_time,
                all_of(names(extra_cols))) %>%
        filter(!is.na(PMID)) %>%
        filter(!is.na(control_mean) & !is.na(test_mean)) 

    return(res)
}

get_data("FST") %>% 
bind_rows(get_data("SPT", 
    c("water_depriv_h" = "Duration.of.prior.water.deprivation..h.",
    "food_depriv_h" = "Duration.of.prior.food.deprivation..h."))) %>%
bind_rows(get_data("EPM")) %>%
bind_rows(get_data("OFT")) %>%
write.csv("CS_data.csv", row.names = FALSE)
