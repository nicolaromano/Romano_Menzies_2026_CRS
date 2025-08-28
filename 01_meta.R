library(metafor)
library(dplyr)

CS_data <- read.csv("CS_data.csv")

get_I2 <- function(model) {
    # See https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
    W <- diag(1 / model$vi)
    X <- model.matrix(model)
    P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
    I2 <- 100 * sum(model$sigma2) / (sum(model$sigma2) +
        (model$k - model$p) / sum(diag(P)))
    return(I2)
}

do_meta <- function(
    test, save_pdf = FALSE, effects_limits = NA) {
    # Calculate effect sizes for the specified test
    # We use Hedges' g for standardized mean difference
    effects <- escalc(
        data = CS_data %>%
            filter(Outcome == test),
        measure = "SMD",
        m1i = control_mean,
        m2i = test_mean,
        n1i = control_n,
        n2i = test_n,
        sd1i = control_sd,
        sd2i = test_sd,
        append = TRUE, # Return the data along with the effect sizes
        slab = PubID
        #slab = paste(PubID, ifelse(is.na(sex), "?", sex), sep = " - ")
    )

    # Mixed-effects model
    mc <- rma.mv(
        data = effects,
        yi = yi,
        V = vi,
        mods = total_time ~ 1, # We include total_time as a moderator
        # We use PubID as a random effect. This is important as some
        # studies have multiple observations (e.g. M/F or different lengths
        # of CRS).
        random = ~ 1 | PubID,
        method = "REML", # Restricted Maximum Likelihood
        test = "knha" # Knapp-Hartung adjustment for small sample sizes
    )

    print(summary(mc))

    if (!is.null(out_pdf)) {
        pdf(out_pdf, width = width, height = height)
    }

    forest(mc,
        order = total_days, #PubID,
        header = c("Study", "Hedges' g [95% CI]"),
        shade = TRUE,
        ilab = cbind(ifelse(is.na(effects$total_days), "?", effects$total_days),
            ifelse(is.na(effects$sex), "?", effects$sex)),
        ilab.lab = c("CRS days", "Sex"),
        ilab.xpos = c(effects_limits[1] - 2, effects_limits[1]),
        xlim = c(effects_limits[1] - 10, effects_limits[2] + 8),
        at = seq(effects_limits[1], effects_limits[2], 5),
        mlab = paste("Random-effects model", test, sep = " - ")
    )

    mtext(bquote(I^2 == .(format(get_I2(mc), digits = 3))),
        side = 1, line = 0, at = (effects_limits[1] - 9.5),
        cex = 1, adj = 0
    )

    summary(mc)

    if (!is.null(out_pdf)) {
        dev.off()
    }

    return(mc)
}

spt_mc <- do_meta("SPT", effects_limits = c(-5, 15))
# mc_FST <- do_meta("FST", out_pdf = "FST_forest.pdf")
# metareg(mc_FST, ~total_time)
# funnel(mc_FST)
# mc_SPT <- do_meta("SPT", out_pdf = "SPT_forest.pdf", 10, 14)
# metareg(mc_SPT, ~ water_depriv_h + food_depriv_h + total_time)

# mc_EPM <- do_meta("EPM", out_pdf = "EPM_forest.pdf")

# mc_OFT <- do_meta("OFT", out_pdf = "OFT_forest.pdf")
