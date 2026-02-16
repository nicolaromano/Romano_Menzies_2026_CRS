library(metafor)
library(dplyr)
library(ggplot2)

out_dir <- "plots"
CRS_data <- read.csv("CRS_data.csv")
predicted_effects <- read.csv("model_estimates.csv")
rownames(predicted_effects) <- predicted_effects$Outcome

run_model <- function(CRS_data, test) {
    effect_sizes <- escalc(
        data = CRS_data,
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

    median_CRS_total_time_h <- predicted_effects[test, "median_CRS_time_h"]

    if (test != "SPT") {
        pred <- predict(model, newmods = median_CRS_total_time_h)
    } else {
        pred <- predict(model, newmods = c(median_CRS_total_time_h, 0))
    }

    return(list(model = model, prediction = pred))
}


# metafor::leave1out does not support rma.mv objects
# See discussion here https://stats.stackexchange.com/questions/155693/
# It makes sense to remove one study at a time (which might have multiple
# effects), not one effect at a time to assess influence yet keep the
# multilevel structure of the model.
leave_one_out <- function(CRS_data, test) {
    CRS_data_filtered <- CRS_data %>%
        filter(Outcome == test) %>%
        filter(!is.na(control_n) & !is.na(test_n)) %>% # Must have sample sizes
        filter(!is.na(control_sd) & !is.na(test_sd)) %>% # Must have SDs
        filter(!is.na(total_time_h)) # Must have CRS duration

    if (test == "SPT") {
        CRS_data_filtered <- CRS_data_filtered %>%
            filter(!is.na(water_depriv_h)) # Must have water deprivation hours for SPT
    }

    clusters <- CRS_data_filtered$PubID %>% unique()

    res <- data.frame(
        PubID = clusters,
        beta_CRS_length = NA,
        SE_CRS_length = NA,
        pval_CRS_length = NA,
        beta_water_deprivation = NA,
        SE_water_deprivation = NA,
        pval_water_deprivation = NA,
        predicted_effect = NA,
        CI_LB = NA,
        CI_UB = NA
    )

    for (i in seq_along(clusters)) {
        CRS_data_filtered_cluster <- CRS_data_filtered %>%
            filter(PubID != clusters[i])

        model_cluster <- run_model(CRS_data_filtered_cluster, test)

        model <- model_cluster$model
        pred <- model_cluster$prediction

        res$beta_CRS_length[i] <- model$beta[2]
        res$SE_CRS_length[i] <- model$se[2]
        res$pval_CRS_length[i] <- model$pval[2]
        if (test == "SPT") {
            res$beta_water_deprivation[i] <- model$beta[3]
            res$SE_water_deprivation[i] <- model$se[3]
            res$pval_water_deprivation[i] <- model$pval[3]
        }
        res$predicted_effect[i] <- pred$pred
        res$CI_LB[i] <- pred$ci.lb
        res$CI_UB[i] <- pred$ci.ub
    }

    return(res)
}

SYRCLE_sensitivity <- function(CRS_data, test) {
    CRS_data_filtered <- CRS_data %>%
        filter(Outcome == test) %>%
        filter(!is.na(control_n) & !is.na(test_n)) %>% # Must have sample sizes
        filter(!is.na(control_sd) & !is.na(test_sd)) %>% # Must have SDs
        filter(!is.na(total_time_h)) # Must have CRS duration

    if (test == "SPT") {
        CRS_data_filtered <- CRS_data_filtered %>%
            filter(!is.na(water_depriv_h)) # Must have water deprivation hours for SPT
    }

    CRS_data_filtered <- CRS_data_filtered %>%
        filter(SYRCLE_Q8 != "No")

    model_res <- run_model(CRS_data_filtered, test)

    return(data.frame(
        PubID = "SYRCLE",
        beta_CRS_length = model_res$model$beta[2],
        SE_CRS_length = model_res$model$se[2],
        pval_CRS_length = model_res$model$pval[2],
        beta_water_deprivation = if (test == "SPT") {
            model_res$model$beta[3]
        } else {
            NA
        },
        SE_water_deprivation = if (test == "SPT") {
            model_res$model$se[3]
        } else {
            NA
        },
        pval_water_deprivation = if (test == "SPT") {
            model_res$model$pval[3]
        } else {
            NA
        },
        predicted_effect = model_res$pred$pred,
        CI_LB = model_res$pred$ci.lb,
        CI_UB = model_res$pred$ci.ub
    ))
}

plot_sensitivity_analysis <- function(outcome, type = "effect", xlim, outfile = NULL) {
    stopifnot(type %in% c("effect", "beta"))

    loo <- leave_one_out(CRS_data, outcome)

    loo$perc_diff <- (predicted_effects[outcome, "predicted_effect"] - loo$predicted_effect) /
        (predicted_effects[outcome, "predicted_effect"]) * 100

    # Order studies by effect
    if (type == "effect") {
        loo <- loo %>%
            arrange(predicted_effect)
    } else {
        loo <- loo %>%
            arrange(beta_CRS_length)
    }

    SYRCLE <- SYRCLE_sensitivity(CRS_data, outcome)
    SYRCLE$perc_diff <- (predicted_effects[outcome, "predicted_effect"] - SYRCLE$predicted_effect) /
        (predicted_effects[outcome, "predicted_effect"]) * 100

    empty_row <- SYRCLE[1, ]
    empty_row[,] <- NA
    empty_row$PubID <- " "
    
    sensit <- rbind(
        SYRCLE,
        empty_row,
        loo
    ) %>%
        mutate(PubID = factor(PubID, levels = c("SYRCLE", " ", loo$PubID)))

    if (!is.null(outfile)) {
        cairo_pdf(outfile, width = 8, height = 8)
    }

    pl <- NULL

    if (type == "effect") {
        pl <- ggplot(sensit, aes(x = predicted_effect, y = PubID)) +
            geom_point(
                size = 3,
                col = rep(c("darkred", "black"), c(2, nrow(loo))),
                pch = rep(c(17, 20), c(2, nrow(loo)))
            ) +
            geom_errorbar(aes(xmin = CI_LB, xmax = CI_UB), width = 0.2) +
            geom_vline(
                xintercept = predicted_effects[outcome, "predicted_effect"],
                linetype = "dashed"
            ) +
            geom_vline(xintercept = 0, linetype = "dotted") +
            geom_rect(
                data = predicted_effects %>% filter(Outcome == outcome),
                aes(xmin = CI_LB, xmax = CI_UB, ymin = -Inf, ymax = Inf),
                fill = "gray", alpha = 0.2, inherit.aes = FALSE
            ) +
            xlim(xlim) +
            ylab("Study removed") +
            xlab("Effect size (Hedges' g)") +
            ggtitle(paste("Sensitivity analysis - ", outcome, "- estimated effect")) +
            theme_minimal(base_size = 14)
    } else { # beta
        pl <- ggplot(sensit, aes(x = beta_CRS_length, y = PubID)) +
            geom_point(
                size = 3,
                col = rep(c("darkred", "black"), c(2, nrow(loo))),
                pch = rep(c(17, 20), c(2, nrow(loo)))
            ) +
            geom_errorbar(aes(
                xmin = beta_CRS_length - 1.96 * SE_CRS_length,
                xmax = beta_CRS_length + 1.96 * SE_CRS_length
            ), width = 0.2) +
            geom_vline(
                xintercept = predicted_effects[outcome, "beta_CRS_length"],
                linetype = "dashed"
            ) +
            geom_rect(
                data = predicted_effects %>% filter(Outcome == outcome),
                aes(
                    xmin = beta_CRS_length - 1.96 * SE_CRS_length,
                    xmax = beta_CRS_length + 1.96 * SE_CRS_length,
                    ymin = -Inf, ymax = Inf
                ),
                fill = "gray", alpha = 0.2, inherit.aes = FALSE
            ) +
            geom_vline(xintercept = 0, linetype = "dotted") +
            xlim(xlim) +
            ylab("Study removed") +
            # xlab(expression(\beta_{"CRS length"})) +
            ggtitle(
                bquote("Sensitivity analysis - " ~ .(outcome) ~
                    " - " ~ beta[CRS ~ length])
            ) +
            theme_minimal(base_size = 14)
    }

    print(pl)

    if (!is.null(outfile)) {
        dev.off()
    }

    return(pl)
}

plot_sensitivity_analysis("FST", "effect", xlim = c(0, 4), paste0(out_dir, "/FST_sensitivity_effect.pdf"))
plot_sensitivity_analysis("SPT", "effect", xlim = c(-5, 1), paste0(out_dir, "/SPT_sensitivity_effect.pdf"))
plot_sensitivity_analysis("EPM", "effect", xlim = c(-3, 1), paste0(out_dir, "/EPM_sensitivity_effect.pdf"))
plot_sensitivity_analysis("OFT", "effect", xlim = c(-4, 4), paste0(out_dir, "/OFT_sensitivity_effect.pdf"))

plot_sensitivity_analysis("FST", "beta", xlim = c(-0.04, 0.04), paste0(out_dir, "/FST_sensitivity_beta.pdf"))
plot_sensitivity_analysis("SPT", "beta", xlim = c(-0.03, 0.01), paste0(out_dir, "/SPT_sensitivity_beta.pdf"))
plot_sensitivity_analysis("EPM", "beta", xlim = c(-0.03, 0.02), paste0(out_dir, "/EPM_sensitivity_beta.pdf"))
plot_sensitivity_analysis("OFT", "beta", xlim = c(-0.05, 0.05), paste0(out_dir, "/OFT_sensitivity_beta.pdf"))
