library(meta)
library(dplyr)

CS_data <- read.csv("CS_data.csv")

do_meta <- function(test, out_pdf = NULL, width = 10, height = 10) {
    mc <- metacont(
        data = CS_data %>%
            filter(Outcome == test),
        n.c = control_n,
        mean.c = control_mean,
        sd.c = control_sd,
        n.e = test_n,
        mean.e = test_mean,
        sd.e = test_sd,
        method.smd = "Hedges", # Hedges' g for standardized mean difference
        studlab = PubID, # Study labels
        common = FALSE,
        random = TRUE, # Random effects model
        # Hartung-Knapp adjustment to improve precision of CI
        method.random.ci = "HK",
        # Method for heterogeneity estimation
        # "REML" (= Restricted Maximum Likelihood) is
        # recommended for continuous outcomes
        method.tau = "REML",
        prediction = TRUE,
        title = test
    )

    print(summary(mc))

    if (!is.null(out_pdf)) {
        pdf(out_pdf, width = width, height = height)
    }

    meta::forest(mc,
        sortvar = PubID, # Sort by publication ID
        xlab = "Standardized Mean Difference (Hedges' g)",
        fontsize = 8,
        col.square = "gray", # Color for squares
        col.diamond = "darkgray", # Color for diamond
        col.study = "black", # Color for study labels
        col.inside = "black",
        leftcols = c("studlab"),
        title = (test),
        plotwidth = "7cm"
    )

    if (!is.null(out_pdf)) {
        dev.off()
    }

    return(mc)
}

mc_FST <- do_meta("FST", out_pdf = "FST_forest.pdf")
metareg(mc_FST, ~total_time)
funnel(mc_FST)
mc_SPT <- do_meta("SPT", out_pdf = "SPT_forest.pdf", 10, 14)
metareg(mc_SPT, ~ water_depriv_h + food_depriv_h + total_time)

mc_EPM <- do_meta("EPM", out_pdf = "EPM_forest.pdf")

mc_OFT <- do_meta("OFT", out_pdf = "OFT_forest.pdf")
