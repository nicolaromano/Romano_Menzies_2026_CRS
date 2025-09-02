library(metafor)
library(dplyr)
library(ggplot2)

CS_data <- read.csv("CS_data.csv")

CS_data %>% 
    group_by(Outcome) %>%
    summarise(total_days_not_specified = sum(is.na(total_days)),
              total_days_min = min(total_days, na.rm = TRUE),
              total_days_max = max(total_days, na.rm = TRUE),
              n_studies = n_distinct(PubID),
              n_effects = n())

cairo_pdf("CRS_duration.pdf",
    width = 10, height = 10,
    family = "Noto Sans"
)

CS_data %>%
    filter(!is.na(total_days)) %>%
    group_by(Outcome) %>%
    select(PMID, total_days, Outcome) %>%
    distinct() %>%
    ggplot(aes(y = total_days, x = Outcome)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, height = 0, size = 2) +
    theme_minimal() +
    ylab("CRS duration (days)") +
    xlab("Behavioral test") +
    theme(text = element_text(family = "Noto Sans"),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14))

dev.off()

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
    effect_sizes <- escalc(
        data = CS_data %>%
            filter(Outcome == test),
        measure = "SMD",
        m1i = test_mean,
        m2i = control_mean,
        n1i = test_n,
        n2i = control_n,
        sd1i = test_sd,
        sd2i = control_sd,
        append = TRUE, # Return the data along with the effect sizes
        slab = PubID
        # slab = paste(PubID, ifelse(is.na(sex), "?", sex), sep = " - ")
    )

    if (is.na(effects_limits[1])) {
        effects_limits <- c(
            floor(min(effect_sizes$yi, na.rm = TRUE)) - 1,
            ceiling(max(effect_sizes$yi, na.rm = TRUE)) + 1
        )
    }

    # Mixed-effects model
    mc <- rma.mv(
        data = effect_sizes,
        yi = yi,
        V = vi,
        mods = ~ total_time, # We include total_time as a moderator
        # We use PubID as a random effect. This is important as some
        # studies have multiple observations (e.g. M/F or different lengths
        # of CRS).
        random = ~ 1 | PubID,
        method = "REML", # Restricted Maximum Likelihood
        test = "knha" # Knapp-Hartung adjustment for small sample sizes
    )

    print(summary(mc))

    if (save_pdf) {
        out_pdf <- paste0(test, "_forest.pdf")
        cairo_pdf(out_pdf,
            width = 10, height = 2 * nrow(effect_sizes) / 5,
            family = "Noto Sans"
        )
    }

    par(family = "Noto Sans")

    if (test == "SPT") {
        ilabs <- cbind(
            ifelse(is.na(effect_sizes$total_days), "?", effect_sizes$total_days),
            ifelse(is.na(effect_sizes$water_depriv_h), "?", effect_sizes$water_depriv_h),
            ifelse(is.na(effect_sizes$sex), "?", effect_sizes$sex)
        )
        ilabs_labels <- c("CRS\ndays", "Water\ndepriv (h)", "Sex")
        ilabs_xpos <- c(
            effects_limits[1] - 4, effects_limits[1] - 2,
            effects_limits[1]
        )
    } else {
        ilabs <- cbind(
            ifelse(is.na(effect_sizes$total_days), "?", effect_sizes$total_days),
            ifelse(is.na(effect_sizes$sex), "?", effect_sizes$sex)
        )
        ilabs_labels <- c("CRS\ndays", "Sex")
        ilabs_xpos <- c(
            effects_limits[1] - 3, effects_limits[1]
        )
    }
    forest(mc,
        order = total_days, # PubID,
        header = c("Study", "Hedges' g [95% CI]"),
        shade = TRUE,
        ilab = ilabs,
        ilab.lab = ilabs_labels,
        ilab.xpos = ilabs_xpos,
        xlim = c(effects_limits[1] - 10, effects_limits[2] + 8),
        ylim = c(-3, nrow(effect_sizes) + 3),
        at = seq(effects_limits[1], effects_limits[2], 5),
        mlab = paste("Random-effects model", test, sep = " - "),
        addfit = FALSE
    )

    median_CRS_days <- median(effect_sizes$total_days, na.rm = TRUE)
    pred <- predict(mc, newmods = median_CRS_days)
    estimate <- pred$pred
    se  <- pred$se
    df  <- mc$ddf[1]

    tval <- estimate / se
    pval <- 2 * (1 - pt(abs(tval), df))

    print(paste0("Model estimate for ", test, 
        " at median CRS days (", median_CRS_days, ")"))
    print(paste0("Hedges' g = ", round(estimate, 3),
                 " (SE = ", round(se, 3),
                 ", df = ", round(df, 1),
                 ", t = ", round(tval, 3),
                 ", p = ", signif(pval, 3), ")"))

    addpoly(
    pred,
    row = -1,
    mlab = "Model estimate (at median CRS days)"
    )

    mtext(bquote(I^2 == .(format(get_I2(mc), digits = 3))),
        side = 1, line = 0, at = (effects_limits[1] - 9.5),
        cex = 1, adj = 0
    )

    summary(mc)

    if (save_pdf) {
        dev.off()
    }

    if (save_pdf) {
        out_pdf <- paste0(test, "_funnel.pdf")
        cairo_pdf(out_pdf,
            width = 7, height = 7,
            family = "Noto Sans"
        )
    }

    # Note: metafor::regtest does not support rma.mv objects
    # Egger's test is a simple linear regression of the effect sizes
    # on their standard errors, weighted by the inverse of the
    # variance of the effect sizes.
    egger <- lm(yi ~ sei, weights = 1/vi, 
        data = transform(effect_sizes, sei = sqrt(vi)))
    print(paste("Egger's test for funnel plot asymmetry - ", test))
    print(summary(egger))

    # Funnel plot
    funnel(mc,
        xlab = "Hedges' g",
        ylab = "Standard Error",
        xlim = effects_limits + c(-5, 5),
        ylim = c(0, max(sqrt(mc$vi)) + 0.1),
        pch = 21,
        las = 1
    )

    intercept <- round(coef(egger)[1], 2)
    pval <- signif(summary(egger)$coefficients[1,4], 2)
    text(x = min(effect_sizes$yi), 
     y = 0, 
     labels = paste0("Egger's intercept = ", intercept, 
                     ", p = ", pval), 
     pos = 4, cex=0.8)

    if (save_pdf) {
        dev.off()
    }

    return(mc)
}

spt_mc <- do_meta("SPT", effects_limits = c(-5, 15), save_pdf = TRUE)
fst_mc <- do_meta("FST", effects_limits = c(-5, 10), save_pdf = TRUE)
epm_mc <- do_meta("EPM", effects_limits = c(-5, 10), save_pdf = TRUE)
oft_mc <- do_meta("OFT", effects_limits = c(-15, 15), save_pdf = TRUE)

# funnel(mc_FST)
# mc_SPT <- do_meta("SPT", out_pdf = "SPT_forest.pdf", 10, 14)
# metareg(mc_SPT, ~ water_depriv_h + food_depriv_h + total_time)

# mc_EPM <- do_meta("EPM", out_pdf = "EPM_forest.pdf")

# mc_OFT <- do_meta("OFT", out_pdf = "OFT_forest.pdf")
