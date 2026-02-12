library(metafor)
library(dplyr)
library(ggplot2)

out_dir <- "plots"
CRS_data <- read.csv("CRS_data.csv")

if (!dir.exists(out_dir)) {
    dir.create(out_dir)
}

CRS_data_summary <- CRS_data %>%
    group_by(Outcome) %>%
    summarise(
        total_days_min = min(total_days, na.rm = TRUE),
        total_days_max = max(total_days, na.rm = TRUE),
        total_hours_min = min(total_time_h, na.rm = TRUE),
        total_hours_max = max(total_time_h, na.rm = TRUE),
        n_studies = n_distinct(PubID),
        # Number of effects is number of rows in the data for that outcome
        # We need to remove rows with NA in control_n or test_n as they
        # cannot be used for meta-analysis
        n_excluded = sum(is.na(control_n) | is.na(test_n) |
            is.na(control_sd) | is.na(test_sd) |
            is.na(total_time_h) |
            (Outcome == "SPT" & is.na(water_depriv_h))),
        n_effects = n() - n_excluded
    )
print(CRS_data_summary)

CRS_data %>%
    group_by(Outcome) %>%
    summarise(
        n_excl = sum(is.na(control_n) | is.na(test_n) |
            is.na(control_sd) | is.na(test_sd) | is.na(total_time_h) | is.na(water_depriv_h)),
        n_studies_excluded = n_distinct(PubID[is.na(control_sd) | is.na(test_sd) | is.na(control_n) | is.na(test_n) |
            is.na(total_time_h) | (Outcome == "SPT" & is.na(water_depriv_h))])
    )

png(paste0(out_dir, "/CRS_time_vs_days.png"),
    width = 10, height = 8, units = "in", res = 300
)

CRS_data %>%
    select(PMID, total_time_h, total_days, Outcome) %>%
    distinct(PMID, .keep_all = TRUE) %>%
    filter(!is.na(total_time_h) & !is.na(total_days)) %>%
    ggplot(aes(x = total_time_h, y = total_days)) +
    geom_point(size = 2) +
    stat_smooth(method = "lm", se = TRUE, color = "#00b7ff", alpha = 0.2) +
    xlab("Total CRS time (hours)") +
    ylab("Days of CRS") +
    theme_minimal() +
    theme(
        axis.text = element_text(family = "Noto Sans", size = 14),
        axis.title = element_text(family = "Noto Sans", size = 16)
    )

dev.off()

CRS_data %>%
    select(PMID, total_time_h, total_days, Outcome) %>%
    distinct(PMID, .keep_all = TRUE) %>%
    filter(!is.na(total_time_h) & !is.na(total_days)) %>%
    lm(total_days ~ total_time_h, data = .) %>%
    summary()

get_I2 <- function(model) {
    # See https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
    W <- diag(1 / model$vi)
    X <- model.matrix(model)
    P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
    I2 <- 100 * sum(model$sigma2) / (sum(model$sigma2) +
        (model$k - model$p) / sum(diag(P)))
    return(I2)
}

add_diamond <- function(model, test, moderators_val, row = -1) {
    pred <- predict(model, newmods = moderators_val)
    estimate <- pred$pred
    se <- pred$se
    ci_lb <- pred$ci.lb
    ci_ub <- pred$ci.ub
    df <- model$ddf[1]

    tval <- estimate / se
    pval <- 2 * (1 - pt(abs(tval), df))

    print(moderators_val)

    print(paste0(
        "Hedges' g = ", round(estimate, 3),
        " (SE = ", round(se, 3),
        ", 95% CI [", round(ci_lb, 3), ", ", round(ci_ub, 3), "]",
        ", df = ", round(df, 1),
        ", t = ", round(tval, 3),
        ", p = ", signif(pval, 3), ")"
    ))

    if (test == "SPT") {
        diamond_label <- paste(
            "Model estimate (at", moderators_val["total_time_h"],
            "h CRS,", moderators_val["water_depriv_h"], "h H2O depr.)"
        )
    } else {
        diamond_label <- paste("Model estimate (at", moderators_val["total_time_h"], "h CRS)")
    }

    addpoly(
        pred,
        row = row,
        mlab = diamond_label
    )
}

CRS_data %>%
    group_by(Outcome) %>%
    summarise(
        median_duration_h = median(total_time_h, na.rm = TRUE)
    )

do_meta <- function(
  test, effects_limits = NA,
  save_pdf = FALSE, pdf_width, pdf_height,
  excluded_studies = NULL
) {
    CRS_data_filtered <- CRS_data %>%
        filter(Outcome == test) %>%
        filter(!(PMID %in% excluded_studies)) %>%
        filter(!is.na(control_n) & !is.na(test_n)) %>% # Must have sample sizes
        filter(!is.na(control_sd) & !is.na(test_sd)) %>% # Must have SDs
        filter(!is.na(total_time_h)) # Must have CRS duration

    if (test == "SPT") {
        CRS_data_filtered <- CRS_data_filtered %>%
            filter(!is.na(water_depriv_h)) # Must have water deprivation hours for SPT
    }

    # Calculate effect sizes for the specified test
    # We use Hedges' g for standardized mean difference
    effect_sizes <- escalc(
        data = CRS_data_filtered,
        measure = "SMD",
        m1i = test_mean,
        m2i = control_mean,
        n1i = test_n,
        n2i = control_n,
        sd1i = test_sd,
        sd2i = control_sd,
        append = TRUE, # Return the data along with the effect sizes
        slab = PubID
    )

    if (is.na(effects_limits[1])) {
        effects_limits <- c(
            floor(min(effect_sizes$yi, na.rm = TRUE)) - 1,
            ceiling(max(effect_sizes$yi, na.rm = TRUE)) + 1
        )
    }

    # Mixed-effects model
    model <- NULL
    if (test != "SPT") {
        model <- rma.mv(
            data = effect_sizes,
            yi = yi,
            V = vi,
            mods = ~total_time_h, # We include total_time as a moderator
            # We use PubID as a random effect. This is important as some
            # studies have multiple observations (e.g. M/F or different lengths
            # of CRS).
            random = ~ 1 | PubID,
            method = "REML", # Restricted Maximum Likelihood
            test = "knha" # Knapp-Hartung adjustment for small sample sizes
        )
    } else {
        model <- rma.mv(
            data = effect_sizes,
            yi = yi,
            V = vi,
            # For SPT, we also include water deprivation as a moderator
            mods = ~ total_time_h + water_depriv_h,
            random = ~ 1 | PubID,
            method = "REML",
            test = "knha"
        )
    }

    print(summary(model))

    if (save_pdf) {
        out_pdf <- paste0(out_dir, "/", test, "_forest.pdf")
        cairo_pdf(out_pdf,
            width = pdf_width, height = pdf_height,
            family = "Noto Sans"
        )
    }

    par(family = "Noto Sans")

    if (test == "SPT") {
        ilabs <- cbind(
            effect_sizes$total_time_h,
            effect_sizes$water_depriv_h,
            ifelse(is.na(effect_sizes$sex), "?", effect_sizes$sex),
            round(weights(model), 1)
        )
        ilabs_labels <- c("Tot.\nCRS (h)", "Water\ndepr. (h)", "Sex", "Study\nweight %")
        ilabs_xpos <- c(
            effects_limits[1] - 4, effects_limits[1] - 2,
            effects_limits[1], effects_limits[2] + 1
        )
    } else {
        ilabs <- cbind(
            ifelse(is.na(effect_sizes$total_time_h), "?", effect_sizes$total_time_h),
            ifelse(is.na(effect_sizes$sex), "?", effect_sizes$sex),
            round(weights(model), 1)
        )
        ilabs_labels <- c("Tot.\nCRS (h)", "Sex", "Study\nweight %")
        ilabs_xpos <- c(
            effects_limits[1] - 3, effects_limits[1],
            effects_limits[2] + 2
        )
    }

    if (test != "SPT") {
        xlims <- c(effects_limits[1] - 10, effects_limits[2] + 12)
        ylims <- c(-3, model$k + 3) # Number of effects + space for header and model summary
    } else {
        xlims <- c(effects_limits[1] - 15, effects_limits[2] + 14)
        ylims <- c(-4, model$k + 4) # Extra space for second diamond
    }

    forest(model,
        order = paste(
            sprintf("%03d", total_time_h),
            sprintf("%03d", water_depriv_h)
        ),
        header = c("Study", "Hedges' g [95% CI]"),
        shade = TRUE,
        ilab = ilabs,
        ilab.lab = ilabs_labels,
        ilab.xpos = ilabs_xpos,
        xlim = xlims,
        ylim = ylims,
        at = seq(effects_limits[1], effects_limits[2], 5),
        mlab = paste("Random-effects model", test, sep = " - "),
        addfit = FALSE,
        efac = 0.75
    )

    median_CRS_total_time_h <- median(effect_sizes$total_time_h, na.rm = TRUE)

    cat("Median CRS total time (hours): ", median_CRS_total_time_h, "\n")

    if (test == "SPT") {
        # For SPT, we set water_depriv_h to 0 (no water deprivation)
        moderators_val <- c(
            total_time_h = median_CRS_total_time_h,
            water_depriv_h = 0
        )
        moderators_val_max <- c(
            total_time_h = median_CRS_total_time_h,
            water_depriv_h = max(effect_sizes$water_depriv_h, na.rm = TRUE)
        )
    } else {
        moderators_val <- c(
            total_time_h = median_CRS_total_time_h
        )
    }

    add_diamond(model, test, moderators_val, row = -1)
    if (test == "SPT") {
        add_diamond(model, test, moderators_val_max, row = -2)
    }

    mtext(bquote(I^2 == .(format(get_I2(model), digits = 3))),
        side = 1, line = -1, at = (effects_limits[1] - 9.5),
        cex = 1, adj = 0
    )

    summary(model)

    if (save_pdf) {
        dev.off()
    }

    if (save_pdf) {
        out_pdf <- paste0(out_dir, "/", test, "_funnel.pdf")
        cairo_pdf(out_pdf,
            width = 10, height = 10,
            family = "Noto Sans"
        )
    }

    # Note: metafor::regtest does not support rma.mv objects
    # Egger's test is a simple linear regression of the effect sizes
    # on their standard errors, weighted by the inverse of the
    # variance of the effect sizes.
    egger <- lm(yi ~ sei,
        weights = 1 / vi,
        data = transform(effect_sizes, sei = sqrt(vi))
    )
    print(paste("Egger's test for funnel plot asymmetry - ", test))
    print(summary(egger))

    # Funnel plot
    # Note from ?metafor::funnel
    # "For (mixed-effects) meta-regression models (i.e., models involving moderators), the plot shows the residuals on the x-axis against their corresponding standard errors".
    # This also means that the reference line is at 0.
    funnel(model,
        xlab = "Residuals",
        ylab = "Standard Error",
        pch = 20,
        cex = 1.2,
        yaxt = "n", # No y axis ticks, we'll add them manually
        hlines = NULL,
        main = paste("Funnel plot -", test)
    )

    # Get default y axis ticks and re-add them manually
    # formatted with 1 decimal place
    yticks <- axTicks(2)
    axis(2,
        at = yticks,
        labels = formatC(yticks, format = "f", digits = 1),
        las = 1
    )

    # Annotate outlier
    if (test == "FST") {
        text(45, 8, "Zhu et al. 2021", cex = 0.7)
    }

    intercept <- round(coef(egger)[1], 2)
    se <- round(summary(egger)$coefficients[1, 2], 2)
    pval <- signif(summary(egger)$coefficients[1, 4], 2)
    tval <- round(summary(egger)$coefficients[1, 3], 2)
    print(paste0(
        "Egger's intercept (SE) = ", intercept,
        " (", se, "), t (", df.residual(egger), ") = ", tval,
        " p = ", pval
    ))

    if (save_pdf) {
        dev.off()
    }

    return(list(model = model, effects = effect_sizes))
}

print("****** FST ******")
fst_res <- do_meta("FST",
    effects_limits = c(-5, 10),
    save_pdf = TRUE, pdf_width = 12, pdf_height = 12
)

fst_res_excl <- do_meta("FST",
    effects_limits = c(-5, 10),
    save_pdf = FALSE, excluded_studies = 33505499
)
print("****** SPT ******")
spt_res <- do_meta("SPT",
    effects_limits = c(-10, 5),
    save_pdf = TRUE, pdf_width = 12, pdf_height = 15
)
print("****** EPM ******")
epm_res <- do_meta("EPM",
    effects_limits = c(-10, 5),
    save_pdf = TRUE, pdf_width = 12, pdf_height = 12
)
print("****** OFT ******")
oft_res <- do_meta("OFT",
    effects_limits = c(-15, 15),
    save_pdf = TRUE, pdf_width = 12, pdf_height = 10
)

# Save predicted effects and betas for sensitivity analyses
fst_pred <- predict(fst_res$model, newmods = median(fst_res$effects$total_time_h, na.rm = TRUE))
spt_pred <- predict(spt_res$model, newmods = c(total_time_h = median(spt_res$effects$total_time_h, na.rm = TRUE), water_depriv_h = 0))
epm_pred <- predict(epm_res$model, newmods = median(epm_res$effects$total_time_h, na.rm = TRUE))
oft_pred <- predict(oft_res$model, newmods = median(oft_res$effects$total_time_h, na.rm = TRUE))


data.frame(
    Outcome = c("FST", "SPT", "EPM", "OFT"),
    beta_CRS_length = c(
        fst_res$model$beta[1],
        spt_res$model$beta[1],
        epm_res$model$beta[1],
        oft_res$model$beta[1]
    ),
    SE_CRS_length = c(
        fst_res$model$se[1],
        spt_res$model$se[1],
        epm_res$model$se[1],
        oft_res$model$se[1]
    ),
    beta_water_deprivation = c(
        NA,
        spt_res$model$beta[2],
        NA,
        NA
    ),
    SE_water_deprivation = c(
        NA,
        spt_res$model$se[2],
        NA,
        NA
    ),
    predicted_effect = c(
        fst_pred$pred,
        spt_pred$pred,
        epm_pred$pred,
        oft_pred$pred
    ),
    CI_LB = c(
        fst_res$model$ci.lb[1],
        spt_res$model$ci.lb[1],
        epm_res$model$ci.lb[1],
        oft_res$model$ci.lb[1]
    ),
    CI_UB = c(
        fst_res$model$ci.ub[1],
        spt_res$model$ci.ub[1],
        epm_res$model$ci.ub[1],
        oft_res$model$ci.ub[1]
    )
) %>%
    write.csv("model_estimates.csv", row.names = FALSE)
