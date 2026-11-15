# Code for the analysis in "Extending multilevel network meta-regression to
# disconnected networks and single-arm studies" (Chapter 3 / paper)

# remotes::install_github("dmphillippo/multinma", ref = "develop")  # run once, see README
library(multinma)
library(tidyverse)
library(ggtext)
library(patchwork)
library(ggraph)
library(loo)

set.seed(12345)

options(mc.cores = parallel::detectCores())

# IPD studies
pso_ipd <- plaque_psoriasis_ipd %>%
  mutate(
    # Variable transformations
    bsa = bsa / 100,
    weight = weight / 10,
    durnpso = durnpso / 10,
    prevsys = as.numeric(prevsys),
    psa = as.numeric(psa),
    age = age / 10,
    pasi_w0 = pasi_w0 / 10,
    # Treatment classes
    trtclass = case_when(trtn == 1 ~ "Placebo",
                         trtn %in% c(2, 3, 5, 6) ~ "IL-17 blocker",
                         trtn == 4 ~ "TNFa blocker",
                         trtn == 7 ~ "IL-12/23 blocker"),
    # Check complete cases for covariates of interest
    is_complete = complete.cases(durnpso, prevsys, bsa, weight, psa)
  ) %>%
  arrange(studyc, trtn)

pso_ipd <- pso_ipd %>%
  mutate(across(c(male, prevsys, psa), as.numeric))

# AgD studies
pso_agd <- plaque_psoriasis_agd %>%
  mutate(
    # Variable transformations
    bsa_mean = bsa_mean / 100,
    bsa_sd = bsa_sd / 100,
    weight_mean = weight_mean / 10,
    weight_sd = weight_sd / 10,
    durnpso_mean = durnpso_mean / 10,
    durnpso_sd = durnpso_sd / 10,
    prevsys = prevsys / 100,
    psa = psa / 100,
    male = male / 100,
    age_mean = age_mean / 10,
    age_sd = age_sd / 10,
    pasi_w0_mean = pasi_w0_mean / 10,
    pasi_w0_sd = pasi_w0_sd / 10,
    # Treatment classes
    trtclass = case_when(trtn == 1 ~ "Placebo",
                         trtn %in% c(2, 3, 5, 6) ~ "IL-17 blocker",
                         trtn == 4 ~ "TNFa blocker",
                         trtn == 7 ~ "IL-12/23 blocker")
  ) %>%
  arrange(studyc, trtn)


# Missing Data
pso_ipd %>%
  group_by(studyc) %>%
  summarise(n_total = n(),
            n_missing = sum(!is_complete),
            pct_missing = mean(!is_complete) * 100)

pso_ipd <- filter(pso_ipd, is_complete)

# Creating the FULL network
pso_net <- combine_network(
  set_ipd(pso_ipd,
          study = studyc,
          trt = trtc,
          r = multi(r0 = 1,
                    PASI75 = pasi75,
                    PASI90 = pasi90,
                    PASI100 = pasi100,
                    type = "ordered", inclusive = TRUE),
          trt_class = trtclass),
  set_agd_arm(pso_agd,
              study = studyc,
              trt = trtc,
              r = multi(r0 = pasi75_n,
                        PASI75 = pasi75_r,
                        PASI90 = pasi90_r,
                        PASI100 = pasi100_r,
                        type = "ordered", inclusive = TRUE),
              trt_class = trtclass)
)

pso_net

class_pal <- c("#D95F02", "#7570B3", "#E7298A", "#E6AB02",
               "#66A61E")

plot(pso_net, weight_nodes = TRUE, weight_edges = TRUE, show_trt_class = TRUE, nudge = 0.1) +
  ggraph::scale_edge_colour_manual("Data",
                                   values = c(AgD = "#113259", IPD = "#55A480"),
                                   guide = guide_legend(override.aes = list(edge_width = 2))) +
  scale_fill_manual("Treatment class",
                    values = class_pal,
                    aesthetics = c("fill", "colour"),
                    guide = guide_legend(override.aes = list(size = 2)))

pso_net <- add_integration(pso_net,
                           durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
                           prevsys = distr(qbern, prob = prevsys),
                           bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
                           weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
                           psa = distr(qbern, prob = psa),
                           age = distr(qgamma, mean = age_mean, sd = age_sd),
                           male = distr(qbern, prob = male),
                           pasi_w0 = distr(qgamma, mean = pasi_w0_mean, sd = pasi_w0_sd),
                           n_int = 64)

## Target population ----
# External Populations
new_agd_means <- tibble::tribble(
  ~study, ~covariate,  ~mean,   ~sd,
  "PsoBest",      "bsa",     24,  20.5,
  "PsoBest",  "durnpso",   18.2,  14.1,
  "PsoBest",  "prevsys",   0.54,    NA,
  "PsoBest",      "psa",  0.207,    NA,
  "PsoBest",   "weight",     85,  19.1,
  "PsoBest", "age", 47.7,14.7,
  "PsoBest", "male", 0.601, NA,
  "PsoBest", "pasi_w0", 14.7, 9.7,
  "PROSPECT",      "bsa",   18.7,  18.4,
  "PROSPECT",  "durnpso",   19.6,  13.5,
  "PROSPECT",  "prevsys", 0.9095,    NA,
  "PROSPECT",      "psa",  0.202,    NA,
  "PROSPECT",   "weight",   87.5,  20.3,
  "PROSPECT", "age", 48.1, 13.7,
  "PROSPECT", "male", 0.624, NA,
  "PROSPECT", "pasi_w0", 17.7, 12.5,
  "Chiricozzi 2019",      "bsa",     23, 16.79,
  "Chiricozzi 2019",  "durnpso",  16.93, 10.82,
  "Chiricozzi 2019",  "prevsys", 0.9061,    NA,
  "Chiricozzi 2019",      "psa", 0.2152,    NA,
  "Chiricozzi 2019",   "weight",   78.3, 15.87,
  "Chiricozzi 2019", "age", 51.9, 14.6,
  "Chiricozzi 2019", "male", 0.682, NA,
  "Chiricozzi 2019", "pasi_w0", 16.65, 8.25
) %>%
  # Tidy up
  pivot_wider(id_cols = study,
              names_from = covariate,
              values_from = c(mean, sd),
              names_glue = "{covariate}_{.value}") %>%
  # Rescale as per analysis
  transmute(study,
            bsa_mean = bsa_mean / 100,
            bsa_sd = bsa_sd / 100,
            weight_mean = weight_mean / 10,
            weight_sd = weight_sd / 10,
            durnpso_mean = durnpso_mean / 10,
            durnpso_sd = durnpso_sd / 10,
            prevsys = prevsys_mean,
            psa = psa_mean,
            age_mean = age_mean / 10,
            age_sd = age_sd / 10,
            male = male_mean,
            pasi_w0_mean = pasi_w0_mean / 10,
            pasi_w0_sd = pasi_w0_sd / 10
  )

new_agd_int <- add_integration(new_agd_means,
                               durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
                               prevsys = distr(qbern, prob = prevsys),
                               bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
                               weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
                               psa = distr(qbern, prob = psa),
                               age = distr(qgamma, mean = age_mean, sd = age_sd),
                               male = distr(qbern, prob = male),
                               pasi_w0 = distr(qgamma, mean = pasi_w0_mean, sd = pasi_w0_sd),
                               cor = pso_net$int_cor,
                               n_int = 64)

new_agd_int <- filter(new_agd_int, study != "PsoBest")

# Shared setup: category k+1 of fitted_ipd = "achieved at least PASI k" (r0 =
# category 1). fitted_ipd is a Stan array[N] vector[K], flattened
# individual-fastest, so array(fitted_flat, dim = c(ndraws, n_ipd, ncat))
# recovers it directly (no transpose needed).
ipd_thresholds <- function(fit) {
  ipd <- fit$network$ipd
  n_ipd <- nrow(ipd)
  ncat <- ncol(ipd$.r)
  fitted_flat <- as.matrix(fit, pars = "fitted_ipd")
  theta <- array(fitted_flat, dim = c(nrow(fitted_flat), n_ipd, ncat))
  lp <- as.array(fit, pars = "lp__")
  log_lik <- as.matrix(fit, pars = "log_lik")[, 1:n_ipd, drop = FALSE]

  at_least <- function(cat_idx) list(
    p = apply(theta[, , cat_idx:ncat, drop = FALSE], c(1, 2), sum),
    y = as.numeric(rowSums(ipd$.r[, cat_idx:ncat, drop = FALSE]) > 0)
  )
  list(n_ipd = n_ipd, ncat = ncat, log_lik = log_lik, lp = lp, at_least = at_least)
}

# LOO cross-validated predictive metrics (RMSE / classification accuracy) per PASI threshold
loo_cutpoint_metrics <- function(fit, cutpoint_names = c("PASI75", "PASI90", "PASI100")) {
  d <- ipd_thresholds(fit)
  stopifnot(d$ncat == length(cutpoint_names) + 1)

  r_eff <- loo::relative_eff(exp(-d$log_lik), chain_id = rep(1:dim(d$lp)[2], each = dim(d$lp)[1]))
  khat <- loo::psis(-d$log_lik, r_eff = r_eff)$diagnostics$pareto_k

  setNames(lapply(seq_along(cutpoint_names), function(k) {
    py <- d$at_least(k + 1L)
    list(
      n = d$n_ipd,
      prevalence = mean(py$y),
      rmse = loo::loo_predictive_metric(py$p, py$y, d$log_lik, r_eff = r_eff, metric = "rmse"),
      acc = loo::loo_predictive_metric(py$p, py$y, d$log_lik, r_eff = r_eff, metric = "acc"),
      pareto_k_max = max(khat),
      pareto_k_prop_over_0.7 = mean(khat > 0.7)
    )
  }), cutpoint_names)
}

# Binary Bayes R2 collapsing the ordinal model to "achieved PASI75 or not" (category >= 2),
# for comparison against the jointly-defined ordinal R2 from bayes_R2()
pasi75_binary_r2 <- function(fit) {
  d <- ipd_thresholds(fit)
  py <- d$at_least(2)
  p <- py$p; y <- py$y; n <- d$n_ipd

  var_fit <- apply(p, 1, var)
  var_res <- apply(p * (1 - p), 1, mean)
  r2_in_sample <- var_fit / (var_fit + var_res)

  psis_obj <- loo::psis(-d$log_lik)
  err_loo <- loo::E_loo(p, psis_obj, log_ratios = -d$log_lik)$value - y

  wts <- matrix(rexp(n * nrow(p)), nrow(p))
  wts <- wts / rowSums(wts)
  wvar <- function(x) (rowSums(sweep(wts, 2, x^2, "*")) - rowSums(sweep(wts, 2, x, "*"))^2) * (n / (n - 1))
  r2_loo <- pmin(pmax(1 - wvar(err_loo) / wvar(y), -1), 1)

  list(in_sample = quantile(r2_in_sample, c(.025, .5, .975)),
       loo = quantile(r2_loo, c(.025, .5, .975)))
}

# Full model

pso_fit_FE_2 <- nma(pso_net,
                  trt_effects = "fixed",
                  link = "probit",
                  regression = ~(durnpso + prevsys + bsa + weight + psa)*.trt + pasi_w0 + age + male,
                  class_interactions = "common",
                  prior_intercept = normal(scale = 10),
                  prior_trt = normal(scale = 10),
                  prior_reg = normal(scale = 10),
                  prior_aux = flat(),
                  QR = TRUE,
                  init_r = 0.2,
                  save_warmup = FALSE)

bayes_R2_FULL <- bayes_R2(pso_fit_FE_2)
bayes_R2_FULL

loo::loo_predictive_metric(pso_fit_FE_2)

(pso_releff_FE <- marginal_effects(pso_fit_FE_2, mtype = "link"))
(pso_pred_FE <- predict(pso_fit_FE_2, type = "response"))

(pso_releff_FE_new <- marginal_effects(pso_fit_FE_2,
                                       mtype = "link",
                                       newdata = new_agd_int,
                                       study = study,
                                       baseline = list(PROSPECT = distr(qbeta, 1156, 1509-1156),
                                                       "Chiricozzi 2019" = distr(qbeta, 243, 330-243)),
                                       baseline_type = "response",
                                       baseline_level = "aggregate"))

(pso_pred_FE_new <- predict(pso_fit_FE_2,
                            type = "response",
                            newdata = new_agd_int,
                            study = study,
                            baseline = list(PROSPECT = distr(qbeta, 1156, 1509-1156),
                                            "Chiricozzi 2019" = distr(qbeta, 243, 330-243)),
                            baseline_type = "response",
                            baseline_level = "aggregate",
                            trt_ref = "SEC_300"))

loo_cutpoint_metrics(pso_fit_FE_2)

rm(pso_fit_FE_2)

pso_fit_FE_UA <- nma(pso_net,
                  trt_effects = "fixed",
                  link = "probit",
                  class_interactions = "common",
                  prior_intercept = normal(scale = 10),
                  prior_trt = normal(scale = 10),
                  prior_reg = normal(scale = 10),
                  prior_aux = flat(),
                  QR = TRUE,
                  init_r = 0.2,
                  save_warmup = FALSE)

bayes_R2_FULL_UA <- bayes_R2(pso_fit_FE_UA)
bayes_R2_FULL_UA

(pso_pred_FE_new_UA_chir <- predict(pso_fit_FE_UA,
                                    type = "response",
                                    baseline = distr(qbeta, 243, 330-243),
                                    baseline_type = "response",
                                    trt_ref = "SEC_300"))

(pso_pred_FE_new_UA_PROS <- predict(pso_fit_FE_UA,
                                    type = "response",
                                    baseline = distr(qbeta, 1156, 1509-1156),
                                    baseline_type = "response",
                                    trt_ref = "SEC_300"))

(pso_releff_FE_UA <- relative_effects(pso_fit_FE_UA))
(pso_pred_FE_UA <- predict(pso_fit_FE_UA, type = "response"))

loo_cutpoint_metrics(pso_fit_FE_UA)

rm(pso_fit_FE_UA)

#Creating the disconnected network ----
pso_agd_filtered <- pso_agd %>%
  filter(
    studyc != "CLEAR",
    !(studyc == "ERASURE" & trtc_long == "Placebo"),
    !(studyc == "FEATURE" & trtc_long == "Placebo"),
    !(studyc == "FIXTURE" & trtc_long == "Secukinumab 150 mg"),
    !(studyc == "FIXTURE" & trtc_long == "Secukinumab 300 mg"),
    !(studyc == "JUNCTURE" & trtc_long == "Placebo")
  )

pso_agd_no_fixture <- pso_agd_filtered %>%
  filter(studyc != "FIXTURE")

## Creating network with single arm study
pso_agd_filtered_single_arm <- pso_agd_filtered %>%
  filter(
    !(studyc == "FIXTURE" & trtc_long == "Placebo"),
  )


# Creating the disconnected network
pso_net_disconnect <- combine_network(
  set_ipd(pso_ipd,
          study = studyc,
          trt = trtc,
          r = multi(r0 = 1,
                    PASI75 = pasi75,
                    PASI90 = pasi90,
                    PASI100 = pasi100,
                    type = "ordered", inclusive = TRUE),
          trt_class = trtclass),
  set_agd_arm(pso_agd_filtered,
              study = studyc,
              trt = trtc,
              r = multi(r0 = pasi75_n,
                        PASI75 = pasi75_r,
                        PASI90 = pasi90_r,
                        PASI100 = pasi100_r,
                        type = "ordered", inclusive = TRUE),
              trt_class = trtclass),
  trt_ref = "PBO")

pso_net_disconnect

plot(pso_net_disconnect, weight_nodes = TRUE, weight_edges = TRUE, show_trt_class = TRUE, nudge = 0.1) +
  ggraph::scale_edge_colour_manual("Data",
                                   values = c(AgD = "#113259", IPD = "#55A480"),
                                   guide = guide_legend(override.aes = list(edge_width = 2))) +
  scale_fill_manual("Treatment class",
                    values = class_pal,
                    aesthetics = c("fill", "colour"),
                    guide = guide_legend(override.aes = list(size = 2)))

pso_net_disconnect <- add_integration(pso_net_disconnect,
                                      durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
                                      prevsys = distr(qbern, prob = prevsys),
                                      bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
                                      weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
                                      psa = distr(qbern, prob = psa),
                                      age = distr(qgamma, mean = age_mean, sd = age_sd),
                                      male = distr(qbern, prob = male),
                                      pasi_w0 = distr(qgamma, mean = pasi_w0_mean, sd = pasi_w0_sd),
                                      n_int = 64)

p <- plot(pso_net_disconnect,
           weight_nodes = TRUE,
           weight_edges = TRUE,
           show_trt_class = TRUE,
           nudge = 0.1) +
     ggraph::scale_edge_colour_manual("Data",
                                      values = c(AgD = "#113259", IPD = "#55A480"),
                                      guide = guide_legend(override.aes = list(edge_width = 2))) +
     scale_fill_manual("Treatment class",
                       values = class_pal,
                       aesthetics = c("fill", "colour"),
                       guide = guide_legend(override.aes = list(size = 2)))

p +
  annotate("label", x = -1, y = 0.5, label = "Subnetwork 1\n(Disconnected)",
           fill = "white", alpha = 0.8) +
  annotate("label", x = 0.5, y = 1.2, label = "Subnetwork 2\n(Reference Network)",
           fill = "white", alpha = 0.8) +
  # Expand the plot limits so the labels don't get cut off
  expand_limits(y = 1.5)

# Create connected network using external control
# External single-arm studies (secukinumab 300 mg), covariates on the analysis scale
pso_agd_PROSPECT <- tibble(
  studyc = "PROSPECT", trtc_long = "Secukinumab 300 mg", trtc = "SEC_300", trtn = 6,
  trtclass = "IL-17 blocker",
  pasi75_r = 1156, pasi75_n = 1509, pasi90_r = 801, pasi90_n = 1509,
  pasi100_r = 337, pasi100_n = 1509, sample_size_w0 = 1509,
  age_mean = 4.81, age_sd = 1.37, male = 0.624,
  pasi_w0_mean = 1.77, pasi_w0_sd = 1.25,
  bsa_mean = 0.187, bsa_sd = 0.184, weight_mean = 8.75, weight_sd = 2.03,
  durnpso_mean = 1.96, durnpso_sd = 1.35, prevsys = 0.9095, psa = 0.202
)

pso_agd_Chiricozzi <- tibble(
  studyc = "Chiricozzi 2019", trtc_long = "Secukinumab 300 mg", trtc = "SEC_300", trtn = 6,
  trtclass = "IL-17 blocker",
  pasi75_r = 243, pasi75_n = 330, pasi90_r = 127, pasi90_n = 330,
  pasi100_r = 71, pasi100_n = 330, sample_size_w0 = 330,
  age_mean = 5.194, age_sd = 1.464, male = 0.6818,
  pasi_w0_mean = 1.665, pasi_w0_sd = 0.825,
  bsa_mean = 0.23, bsa_sd = 0.1679, weight_mean = 7.83, weight_sd = 1.587,
  durnpso_mean = 1.693, durnpso_sd = 1.082, prevsys = 0.9061, psa = 0.2152
)

pso_agd_CON_PROS_single <- bind_rows(pso_agd_filtered, pso_agd_PROSPECT)
pso_agd_CON_chir_single <- bind_rows(pso_agd_filtered, pso_agd_Chiricozzi)

### for comparison table for ALL studies
pso_agd_CON_chir_PROS_single <- rbind(pso_agd_CON_chir_single, pso_agd_CON_PROS_single) %>%
  dplyr::distinct(studyc, trtc, .keep_all = TRUE)

pso_CON_net_chir_PROS_single <- combine_network(
  set_ipd(pso_ipd,
          study = studyc,
          trt = trtc,
          r = multi(r0 = 1,
                    PASI75 = pasi75,
                    PASI90 = pasi90,
                    PASI100 = pasi100,
                    type = "ordered", inclusive = TRUE),
          trt_class = trtclass),
  set_agd_arm(pso_agd_CON_chir_PROS_single,
              study = studyc,
              trt = trtc,
              r = multi(r0 = pasi75_n,
                        PASI75 = pasi75_r,
                        PASI90 = pasi90_r,
                        PASI100 = pasi100_r,
                        type = "ordered", inclusive = TRUE),
              trt_class = trtclass)
)

pso_CON_net_chir_PROS_single



pso_CON_net_chir_PROS_single <- add_integration(pso_CON_net_chir_PROS_single,
                                           durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
                                           prevsys = distr(qbern, prob = prevsys),
                                           bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
                                           weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
                                           psa = distr(qbern, prob = psa),
                                           age = distr(qgamma, mean = age_mean, sd = age_sd),
                                           male = distr(qbern, prob = male),
                                           pasi_w0 = distr(qgamma, mean = pasi_w0_mean, sd = pasi_w0_sd),
                                           n_int = 64)

PROS_compare_ALL <- compare_populations(pso_CON_net_chir_PROS_single,
                                    covariates = c("bsa", "weight", "durnpso", "prevsys", "psa", "age", "male", "pasi_w0"),
                                    method = c("propensity"))

PROS_compare_ALL

## Combining the networks - single arms ----
pso_CON_net_PROS_single <- combine_network(
  set_ipd(pso_ipd,
          study = studyc,
          trt = trtc,
          r = multi(r0 = 1,
                    PASI75 = pasi75,
                    PASI90 = pasi90,
                    PASI100 = pasi100,
                    type = "ordered", inclusive = TRUE),
          trt_class = trtclass),
  set_agd_arm(pso_agd_CON_PROS_single,
              study = studyc,
              trt = trtc,
              r = multi(r0 = pasi75_n,
                        PASI75 = pasi75_r,
                        PASI90 = pasi90_r,
                        PASI100 = pasi100_r,
                        type = "ordered", inclusive = TRUE),
              trt_class = trtclass),
  trt_ref = "PBO"
)

pso_CON_net_PROS_single

pso_CON_net_PROS_single <- add_integration(pso_CON_net_PROS_single,
                                    durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
                                    prevsys = distr(qbern, prob = prevsys),
                                    bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
                                    weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
                                    psa = distr(qbern, prob = psa),
                                    age = distr(qgamma, mean = age_mean, sd = age_sd),
                                    male = distr(qbern, prob = male),
                                    pasi_w0 = distr(qgamma, mean = pasi_w0_mean, sd = pasi_w0_sd),
                                    n_int = 64)

pso_CON_net_chir_single <- combine_network(
  set_ipd(pso_ipd,
          study = studyc,
          trt = trtc,
          r = multi(r0 = 1,
                    PASI75 = pasi75,
                    PASI90 = pasi90,
                    PASI100 = pasi100,
                    type = "ordered", inclusive = TRUE),
          trt_class = trtclass),
  set_agd_arm(pso_agd_CON_chir_single,
              study = studyc,
              trt = trtc,
              r = multi(r0 = pasi75_n,
                        PASI75 = pasi75_r,
                        PASI90 = pasi90_r,
                        PASI100 = pasi100_r,
                        type = "ordered", inclusive = TRUE),
              trt_class = trtclass),
  trt_ref = "PBO"
)

pso_CON_net_chir_single

pso_CON_net_chir_single <- add_integration(pso_CON_net_chir_single,
                                           durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd),
                                           prevsys = distr(qbern, prob = prevsys),
                                           bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd),
                                           weight = distr(qgamma, mean = weight_mean, sd = weight_sd),
                                           psa = distr(qbern, prob = psa),
                                           age = distr(qgamma, mean = age_mean, sd = age_sd),
                                           male = distr(qbern, prob = male),
                                           pasi_w0 = distr(qgamma, mean = pasi_w0_mean, sd = pasi_w0_sd),
                                           n_int = 64)


# Fixed bridge network ----
#Largest ESS%
pso_fit_Stage_two_FIXED <- nma(pso_net_disconnect,
                               trt_effects = "fixed",
                               link = "probit",
                               regression = ~(durnpso + prevsys + bsa + weight + psa)*.trt + pasi_w0 + age + male,
                               class_interactions = "common",
                               connect_baseline = con(type = "fixed",
                                                      studies = c("UNCOVER-3", "ERASURE")),
                               prior_intercept = normal(scale = 10),
                               prior_trt = normal(scale = 10),
                               prior_aux = flat(),
                               QR = TRUE,
                               init_r = 0.2,
                               save_warmup = FALSE)

bayes_R2_fixed <- bayes_R2(pso_fit_Stage_two_FIXED)
loo_cutpoint_metrics(pso_fit_Stage_two_FIXED)

(pso_releff_FE_FIXED <- marginal_effects(pso_fit_Stage_two_FIXED, mtype = "link"))
(pso_pred_FE_FIXED <- predict(pso_fit_Stage_two_FIXED, type = "response"))

(pso_releff_FE_new_FIXED <- marginal_effects(pso_fit_Stage_two_FIXED,
                                             mtype = "link",
                                             newdata = new_agd_int,
                                             study = study,
                                             baseline = list(PROSPECT = distr(qbeta, 1156, 1509-1156),
                                                             "Chiricozzi 2019" = distr(qbeta, 243, 330-243)),
                                             baseline_type = "response",
                                             baseline_level = "aggregate"))

(pso_pred_FE_new_FIXED <- predict(pso_fit_Stage_two_FIXED,
                                  type = "response",
                                  newdata = new_agd_int,
                                  study = study,
                                  baseline = list(PROSPECT = distr(qbeta, 1156, 1509-1156),
                                                  "Chiricozzi 2019" = distr(qbeta, 243, 330-243)),
                                  baseline_type = "response",
                                  baseline_level = "aggregate",
                                  trt_ref = "SEC_300"))

rm(pso_fit_Stage_two_FIXED)

pso_fit_Stage_two_FIXED_UA <- nma(pso_net_disconnect,
                               trt_effects = "fixed",
                               link = "probit",
                               class_interactions = "common",
                               connect_baseline = con(type = "fixed",
                                                      studies = c("UNCOVER-3", "ERASURE")),
                               prior_intercept = normal(scale = 10),
                               prior_trt = normal(scale = 10),
                               prior_reg = normal(scale = 10),
                               prior_aux = flat(),
                               QR = TRUE,
                               init_r = 0.2,
                               save_warmup = FALSE)

loo_cutpoint_metrics(pso_fit_Stage_two_FIXED_UA)
bayes_R2_fixed_UA <- bayes_R2(pso_fit_Stage_two_FIXED_UA)

(pso_releff_FE_FIXED_UA <- relative_effects(pso_fit_Stage_two_FIXED_UA))
(pso_pred_FE_FIXED_UA <- predict(pso_fit_Stage_two_FIXED_UA, type = "response"))

(pso_pred_FE_new_FIXED_UA_PROS <- predict(pso_fit_Stage_two_FIXED_UA,
                                          type = "response",
                                          baseline = distr(qbeta, 1156, 1509-1156),
                                          baseline_type = "response",
                                          trt_ref = "SEC_300"))

(pso_pred_FE_new_FIXED_UA_chir <- predict(pso_fit_Stage_two_FIXED_UA,
                                          type = "response",
                                          baseline = distr(qbeta, 243, 330-243),
                                          baseline_type = "response",
                                          trt_ref = "SEC_300"))

rm(pso_fit_Stage_two_FIXED_UA)

#Smallest ESS%
pso_fit_Stage_two_FIXED_Bad <- nma(pso_net_disconnect,
                               trt_effects = "fixed",
                               link = "probit",
                               regression = ~(durnpso + prevsys + bsa + weight + psa)*.trt + pasi_w0 + age + male,
                               class_interactions = "common",
                               connect_baseline = con(type = "fixed",
                                                      studies = c("JUNCTURE", "IXORA-S")),
                               prior_intercept = normal(scale = 10),
                               prior_trt = normal(scale = 10),
                               prior_reg = normal(scale = 10),
                               prior_aux = flat(),
                               QR = TRUE,
                               init_r = 0.2,
                               save_warmup = FALSE)

bayes_R2_fixed_bad <- bayes_R2(pso_fit_Stage_two_FIXED_Bad)
loo_cutpoint_metrics(pso_fit_Stage_two_FIXED_Bad)

(pso_releff_FE_FIXED_bad <- marginal_effects(pso_fit_Stage_two_FIXED_Bad, mtype = "link"))
(pso_pred_FE_FIXED_bad <- predict(pso_fit_Stage_two_FIXED_Bad, type = "response"))

(pso_releff_FE_new_FIXED_bad <- marginal_effects(pso_fit_Stage_two_FIXED_Bad,
                                                 mtype = "link",
                                                 newdata = new_agd_int,
                                                 study = study,
                                                 baseline = list(PROSPECT = distr(qbeta, 1156, 1509-1156),
                                                                 "Chiricozzi 2019" = distr(qbeta, 243, 330-243)),
                                                 baseline_type = "response",
                                                 baseline_level = "aggregate"))

(pso_pred_FE_new_FIXED_bad <- predict(pso_fit_Stage_two_FIXED_Bad,
                                  type = "response",
                                  newdata = new_agd_int,
                                  study = study,
                                  baseline = list(PROSPECT = distr(qbeta, 1156, 1509-1156),
                                                  "Chiricozzi 2019" = distr(qbeta, 243, 330-243)),
                                  baseline_type = "response",
                                  baseline_level = "aggregate",
                                  trt_ref = "SEC_300"))

pso_fit_Stage_two_FIXED_Bad
rm(pso_fit_Stage_two_FIXED_Bad)

pso_fit_Stage_two_FIXED_Bad_UA <- nma(pso_net_disconnect,
                                  trt_effects = "fixed",
                                  link = "probit",
                                  class_interactions = "common",
                                  connect_baseline = con(type = "fixed",
                                                         studies = c("JUNCTURE", "IXORA-S")),
                                  prior_intercept = normal(scale = 10),
                                  prior_trt = normal(scale = 10),
                                  prior_reg = normal(scale = 10),
                                  prior_aux = flat(),
                                  QR = TRUE,
                                  init_r = 0.2,
                                  save_warmup = FALSE)

bayes_R2_fixed_bad_UA <- bayes_R2(pso_fit_Stage_two_FIXED_Bad_UA)
loo_cutpoint_metrics(pso_fit_Stage_two_FIXED_Bad_UA)

(pso_releff_FE_FIXED_bad_UA <- relative_effects(pso_fit_Stage_two_FIXED_Bad_UA))
(pso_pred_FE_FIXED_bad_UA <- predict(pso_fit_Stage_two_FIXED_Bad_UA, type = "response"))

(pso_pred_FE_new_FIXED_bad_UA_PROS <- predict(pso_fit_Stage_two_FIXED_Bad_UA,
                                           type = "response",
                                           baseline = distr(qbeta, 1156, 1509-1156),
                                           baseline_type = "response",
                                           trt_ref = "SEC_300"))

(pso_pred_FE_new_FIXED_bad_UA_chir <- predict(pso_fit_Stage_two_FIXED_Bad_UA,
                                           type = "response",
                                           baseline = distr(qbeta, 243, 330-243),
                                           baseline_type = "response",
                                           trt_ref = "SEC_300"))

rm(pso_fit_Stage_two_FIXED_Bad_UA)

# Random Bridge ----
## Stage 1: baseline synthesis ----
pso_fit_Stage_one_FE <- baseline_synthesis(pso_net_disconnect,
                                           trt_effects = "fixed",
                                           link = "probit",
                                           regression = ~(durnpso + prevsys + bsa + weight + psa)*.trt + pasi_w0 + age + male,
                                           class_interactions = "common",
                                           prior_intercept = normal(scale = 10),
                                           prior_intercept_sd = half_normal(scale = 5),
                                           prior_trt = normal(scale = 10),
                                           prior_reg = normal(scale = 10),
                                           prior_aux = flat(),
                                           QR = FALSE,
                                           init_r = 0.2,
                                           save_warmup = FALSE)

bl_new <- as.vector(as.array(pso_fit_Stage_one_FE, pars = "baseline_new"))
t_fit <- MASS::fitdistr(bl_new, "t")
t_fit$estimate  # m = location, s = scale, df

# Unadjusted stage one (no regression), used for the _UA_ stage two models
pso_fit_Stage_one_FE_UA <- baseline_synthesis(pso_net_disconnect,
                                              trt_effects = "fixed",
                                              link = "probit",
                                              class_interactions = "common",
                                              prior_intercept = normal(scale = 10),
                                              prior_intercept_sd = half_normal(scale = 5),
                                              prior_trt = normal(scale = 10),
                                              prior_aux = flat(),
                                              QR = FALSE,
                                              init_r = 0.2,
                                              save_warmup = FALSE)

bl_new_UA <- as.vector(as.array(pso_fit_Stage_one_FE_UA, pars = "baseline_new"))
t_fit_UA <- MASS::fitdistr(bl_new_UA, "t")
t_fit_UA$estimate  # m = location, s = scale, df

### 1 Prior ----
pso_fit_Stage_two_Random_1P <- nma(pso_net_disconnect,
                                trt_effects = "fixed",
                                link = "probit",
                                regression = ~(durnpso + prevsys + bsa + weight + psa)*.trt + pasi_w0 + age + male,
                                class_interactions = "common",
                                connect_baseline = list(con(type = "random",
                                                            studies = c("ERASURE"),
                                                            prior_baseline = student_t(t_fit$estimate["m"], t_fit$estimate["s"], t_fit$estimate["df"]))),
                                prior_intercept = normal(scale = 10),
                                prior_trt = normal(scale = 10),
                                prior_reg = normal(scale = 10),
                                prior_aux = flat(),
                                QR = FALSE,
                                init_r = 0.2,
                                save_warmup = FALSE)

bayes_R2_random_1P <- bayes_R2(pso_fit_Stage_two_Random_1P)
loo_cutpoint_metrics(pso_fit_Stage_two_Random_1P)

(pso_releff_FE_RANDOM_1P <- marginal_effects(pso_fit_Stage_two_Random_1P, mtype = "link"))
(pso_pred_FE_RANDOM_1P <- predict(pso_fit_Stage_two_Random_1P, type = "response"))

(pso_releff_FE_new_RANDOM_1P <- marginal_effects(pso_fit_Stage_two_Random_1P,
                                                 mtype = "link",
                                                 newdata = new_agd_int,
                                                 study = study,
                                                 baseline = list(PROSPECT = distr(qbeta, 1156, 1509-1156),
                                                                 "Chiricozzi 2019" = distr(qbeta, 243, 330-243)),
                                                 baseline_type = "response",
                                                 baseline_level = "aggregate"))

(pso_pred_FE_new_RANDOM_1P <- predict(pso_fit_Stage_two_Random_1P,
                                      type = "response",
                                      newdata = new_agd_int,
                                      study = study,
                                      baseline = list(PROSPECT = distr(qbeta, 1156, 1509-1156),
                                                      "Chiricozzi 2019" = distr(qbeta, 243, 330-243)),
                                      baseline_type = "response",
                                      baseline_level = "aggregate",
                                      trt_ref = "SEC_300"))

rm(pso_fit_Stage_two_Random_1P)

pso_fit_Stage_two_Random_UA_1P <- nma(pso_net_disconnect,
                                trt_effects = "fixed",
                                link = "probit",
                                class_interactions = "common",
                                connect_baseline = list(con(type = "random",
                                                            studies = c("ERASURE"),
                                                            prior_baseline = student_t(t_fit_UA$estimate["m"], t_fit_UA$estimate["s"], t_fit_UA$estimate["df"]))),
                                prior_intercept = normal(scale = 10),
                                prior_trt = normal(scale = 10),
                                prior_reg = normal(scale = 10),
                                prior_aux = flat(),
                                QR = FALSE,
                                init_r = 0.2,
                                save_warmup = FALSE)

bayes_R2_random_1P_UA <- bayes_R2(pso_fit_Stage_two_Random_UA_1P)
loo_cutpoint_metrics(pso_fit_Stage_two_Random_UA_1P)

(pso_releff_FE_RANDOM_UA_1P <- relative_effects(pso_fit_Stage_two_Random_UA_1P))
(pso_pred_FE_RANDOM_UA_1P <- predict(pso_fit_Stage_two_Random_UA_1P, type = "response"))

(pso_pred_FE_new_RANDOM_UA_PROS_1P <- predict(pso_fit_Stage_two_Random_UA_1P,
                                              type = "response",
                                              baseline = distr(qbeta, 1156, 1509-1156),
                                              baseline_type = "response",
                                              trt_ref = "SEC_300"))

(pso_pred_FE_new_RANDOM_UA_chir_1P <- predict(pso_fit_Stage_two_Random_UA_1P,
                                              type = "response",
                                              baseline = distr(qbeta, 243, 330-243),
                                              baseline_type = "response",
                                              trt_ref = "SEC_300"))

rm(pso_fit_Stage_two_Random_UA_1P)

### 3 Prior ----
pso_fit_Stage_two_Random_3P <- nma(pso_net_disconnect,
                                   trt_effects = "fixed",
                                   link = "probit",
                                   regression = ~(durnpso + prevsys + bsa + weight + psa)*.trt + pasi_w0 + age + male,
                                   class_interactions = "common",
                                   connect_baseline = list(con(type = "random",
                                                               studies = c("ERASURE", "JUNCTURE", "FEATURE"),
                                                               prior_baseline = student_t(t_fit$estimate["m"], t_fit$estimate["s"], t_fit$estimate["df"]))),
                                   prior_intercept = normal(scale = 10),
                                   prior_trt = normal(scale = 10),
                                   prior_reg = normal(scale = 10),
                                   prior_aux = flat(),
                                   QR = FALSE,
                                   init_r = 0.2,
                                   save_warmup = FALSE)

bayes_R2_random_3P <- bayes_R2(pso_fit_Stage_two_Random_3P)
loo_cutpoint_metrics(pso_fit_Stage_two_Random_3P)

(pso_releff_FE_RANDOM_3P <- marginal_effects(pso_fit_Stage_two_Random_3P, mtype = "link"))
(pso_pred_FE_RANDOM_3P <- predict(pso_fit_Stage_two_Random_3P, type = "response"))

(pso_releff_FE_new_RANDOM_3P <- marginal_effects(pso_fit_Stage_two_Random_3P,
                                                 mtype = "link",
                                                 newdata = new_agd_int,
                                                 study = study,
                                                 baseline = list(PROSPECT = distr(qbeta, 1156, 1509-1156),
                                                                 "Chiricozzi 2019" = distr(qbeta, 243, 330-243)),
                                                 baseline_type = "response",
                                                 baseline_level = "aggregate"))

(pso_pred_FE_new_RANDOM_3P <- predict(pso_fit_Stage_two_Random_3P,
                                      type = "response",
                                      newdata = new_agd_int,
                                      study = study,
                                      baseline = list(PROSPECT = distr(qbeta, 1156, 1509-1156),
                                                      "Chiricozzi 2019" = distr(qbeta, 243, 330-243)),
                                      baseline_type = "response",
                                      baseline_level = "aggregate",
                                      trt_ref = "SEC_300"))

rm(pso_fit_Stage_two_Random_3P)

pso_fit_Stage_two_Random_UA_3P <- nma(pso_net_disconnect,
                                      trt_effects = "fixed",
                                      link = "probit",
                                      class_interactions = "common",
                                      connect_baseline = list(con(type = "random",
                                                                  studies = c("ERASURE", "JUNCTURE", "FEATURE"),
                                                                  prior_baseline = student_t(t_fit_UA$estimate["m"], t_fit_UA$estimate["s"], t_fit_UA$estimate["df"]))),
                                      prior_intercept = normal(scale = 10),
                                      prior_trt = normal(scale = 10),
                                      prior_reg = normal(scale = 10),
                                      prior_aux = flat(),
                                      QR = FALSE,
                                      init_r = 0.2,
                                      save_warmup = FALSE)

bayes_R2_random_3P_UA <- bayes_R2(pso_fit_Stage_two_Random_UA_3P)
loo_cutpoint_metrics(pso_fit_Stage_two_Random_UA_3P)

(pso_releff_FE_RANDOM_UA_3P <- relative_effects(pso_fit_Stage_two_Random_UA_3P))
(pso_pred_FE_RANDOM_UA_3P <- predict(pso_fit_Stage_two_Random_UA_3P, type = "response"))

(pso_pred_FE_new_RANDOM_UA_PROS_3P <- predict(pso_fit_Stage_two_Random_UA_3P,
                                              type = "response",
                                              baseline = distr(qbeta, 1156, 1509-1156),
                                              baseline_type = "response",
                                              trt_ref = "SEC_300"))

(pso_pred_FE_new_RANDOM_UA_chir_3P <- predict(pso_fit_Stage_two_Random_UA_3P,
                                              type = "response",
                                              baseline = distr(qbeta, 243, 330-243),
                                              baseline_type = "response",
                                              trt_ref = "SEC_300"))

rm(pso_fit_Stage_two_Random_UA_3P)

# External control ----

pso_fit_external_PROS <- nma(pso_CON_net_PROS_single,
                             trt_effects = "fixed",
                             link = "probit",
                             regression = ~(durnpso + prevsys + bsa + weight + psa)*.trt + pasi_w0 + age + male,
                             connect_baseline = list(con(type = "fixed",
                                                         studies = c("PROSPECT", "IXORA-S"))),
                             class_interactions = "common",
                             prior_intercept = normal(scale = 10),
                             prior_trt = normal(scale = 10),
                             prior_reg = normal(scale = 10),
                             prior_aux = flat(),
                             QR = TRUE,
                             init_r = 0.2,
                             save_warmup = FALSE)

bayes_R2_external_pros <- bayes_R2(pso_fit_external_PROS)
loo_cutpoint_metrics(pso_fit_external_PROS)

(pso_releff_FE_external_PROS <- marginal_effects(pso_fit_external_PROS, mtype = "link"))
(pso_pred_FE_external_PROS <- predict(pso_fit_external_PROS, type = "response"))

(pso_releff_FE_new_external_PROS <- marginal_effects(pso_fit_external_PROS,
                                                     mtype = "link",
                                                     newdata = new_agd_int,
                                                     study = study,
                                                     baseline = list(PROSPECT = distr(qbeta, 1156, 1509-1156),
                                                                     "Chiricozzi 2019" = distr(qbeta, 243, 330-243)),
                                                     baseline_type = "response",
                                                     baseline_level = "aggregate"))

(pso_pred_FE_new_external_PROS <- predict(pso_fit_external_PROS,
                                          type = "response",
                                          newdata = new_agd_int,
                                          study = study,
                                          baseline = list(PROSPECT = distr(qbeta, 1156, 1509-1156),
                                                          "Chiricozzi 2019" = distr(qbeta, 243, 330-243)),
                                          baseline_type = "response",
                                          baseline_level = "aggregate",
                                          trt_ref = "SEC_300"))

rm(pso_fit_external_PROS)

pso_fit_external_PROS_UA <- nma(pso_CON_net_PROS_single,
                             trt_effects = "fixed",
                             link = "probit",
                             class_interactions = "common",
                             connect_baseline = list(con(type = "fixed",
                                                         studies = c("PROSPECT", "IXORA-S"))),
                             prior_intercept = normal(scale = 10),
                             prior_trt = normal(scale = 10),
                             prior_reg = normal(scale = 10),
                             prior_aux = flat(),
                             QR = TRUE,
                             init_r = 0.2,
                             save_warmup = FALSE)

bayes_R2_external_pros_UA <- bayes_R2(pso_fit_external_PROS_UA)
loo_cutpoint_metrics(pso_fit_external_PROS_UA)

(pso_releff_FE_external_PROS_UA <- relative_effects(pso_fit_external_PROS_UA))
(pso_pred_FE_external_PROS_UA <- predict(pso_fit_external_PROS_UA, type = "response"))

(pso_pred_FE_new_external_PROS_UA_PROS <- predict(pso_fit_external_PROS_UA,
                                                  type = "response",
                                                  baseline = distr(qbeta, 1156, 1509-1156),
                                                  baseline_type = "response",
                                                  trt_ref = "SEC_300"))

(pso_pred_FE_new_external_PROS_UA_chir <- predict(pso_fit_external_PROS_UA,
                                                  type = "response",
                                                  baseline = distr(qbeta, 243, 330-243),
                                                  baseline_type = "response",
                                                  trt_ref = "SEC_300"))

rm(pso_fit_external_PROS_UA)

pso_fit_external_chir <- nma(pso_CON_net_chir_single,
                             trt_effects = "fixed",
                             link = "probit",
                             regression = ~(durnpso + prevsys + bsa + weight + psa)*.trt + pasi_w0 + age + male,
                             connect_baseline = list(con(type = "fixed",
                                                         studies = c("Chiricozzi 2019", "IXORA-S"))),
                             class_interactions = "common",
                             prior_intercept = normal(scale = 10),
                             prior_trt = normal(scale = 10),
                             prior_reg = normal(scale = 10),
                             prior_aux = flat(),
                             QR = TRUE,
                             init_r = 0.2,
                             save_warmup = FALSE)

bayes_R2_external_chir <- bayes_R2(pso_fit_external_chir)
loo_cutpoint_metrics(pso_fit_external_chir)

(pso_releff_FE_external_chir <- marginal_effects(pso_fit_external_chir, mtype = "link"))
(pso_pred_FE_external_chir <- predict(pso_fit_external_chir, type = "response"))

(pso_releff_FE_new_external_chir <- marginal_effects(pso_fit_external_chir,
                                                     mtype = "link",
                                                     newdata = new_agd_int,
                                                     study = study,
                                                     baseline = list(PROSPECT = distr(qbeta, 1156, 1509-1156),
                                                                     "Chiricozzi 2019" = distr(qbeta, 243, 330-243)),
                                                     baseline_type = "response",
                                                     baseline_level = "aggregate"))

(pso_pred_FE_new_external_chir <- predict(pso_fit_external_chir,
                                          type = "response",
                                          newdata = new_agd_int,
                                          study = study,
                                          baseline = list(PROSPECT = distr(qbeta, 1156, 1509-1156),
                                                          "Chiricozzi 2019" = distr(qbeta, 243, 330-243)),
                                          baseline_type = "response",
                                          baseline_level = "aggregate",
                                          trt_ref = "SEC_300"))

rm(pso_fit_external_chir)

pso_fit_external_chir_UA <- nma(pso_CON_net_chir_single,
                                trt_effects = "fixed",
                                link = "probit",
                                connect_baseline = list(con(type = "fixed",
                                                            studies = c("Chiricozzi 2019", "IXORA-S"))),
                                class_interactions = "common",
                                prior_intercept = normal(scale = 10),
                                prior_trt = normal(scale = 10),
                                prior_reg = normal(scale = 10),
                                prior_aux = flat(),
                                QR = TRUE,
                                init_r = 0.2,
                                save_warmup = FALSE)

bayes_R2_external_chir_UA <- bayes_R2(pso_fit_external_chir_UA)
loo_cutpoint_metrics(pso_fit_external_chir_UA)

(pso_releff_FE_external_chir_UA <- relative_effects(pso_fit_external_chir_UA))
(pso_pred_FE_external_chir_UA <- predict(pso_fit_external_chir_UA, type = "response"))

(pso_pred_FE_new_external_chir_UA_PROS <- predict(pso_fit_external_chir_UA,
                                                  type = "response",
                                                  baseline = distr(qbeta, 1156, 1509-1156),
                                                  baseline_type = "response",
                                                  trt_ref = "SEC_300"))

(pso_pred_FE_new_external_chir_UA_chir <- predict(pso_fit_external_chir_UA,
                                                  type = "response",
                                                  baseline = distr(qbeta, 243, 330-243),
                                                  baseline_type = "response",
                                                  trt_ref = "SEC_300"))

rm(pso_fit_external_chir_UA)

# Function to duplicate nma results and add in names of studies. As NMA results are all the same for each study we need to give
# each study these results.
duplicate_and_add_study <- function(df, study_names) {
  duplicated_df <- df
  original_study <- rep(study_names, each = nrow(df))
  duplicated_study <- rep(study_names, each = nrow(duplicated_df))

  df <- df[rep(seq_len(nrow(df)), length(study_names)), ]
  df <- df %>% mutate(.study = original_study)

  duplicated_df <- duplicated_df[rep(seq_len(nrow(duplicated_df)), length(study_names)), ]
  duplicated_df <- duplicated_df %>% mutate(.study = duplicated_study)

  # Combine the original and duplicated datasets
  combined_df <- bind_rows(df, duplicated_df)

  return(combined_df)
}

# Define the study names
study_names_int <- c("IXORA-S", "ERASURE", "FEATURE", "FIXTURE", "JUNCTURE", "UNCOVER-1", "UNCOVER-2", "UNCOVER-3", "Chiricozzi 2019", "PROSPECT","PsoBest")

pasi_categories <- c("PASI75", "PASI90", "PASI100")
add_all_categories <- function(df, categories = pasi_categories) {
  bind_rows(lapply(categories, function(cat) mutate(df, .category = cat)))
}

# Naming Results ----
## Full ----
FULL_inter_probit <- as.data.frame(pso_releff_FE$summary)
FULL_inter_probit_new <- as.data.frame((pso_releff_FE_new))
FULL_inter_probit <- rbind(FULL_inter_probit, FULL_inter_probit_new)
FULL_inter_probit <- FULL_inter_probit %>% mutate(Network = "FULL")
FULL_inter_probit <- FULL_inter_probit %>% mutate(Adjusted = "MLNMR")

FULL_inter_prob <- as.data.frame(pso_pred_FE$summary)
FULL_inter_prob_new <- as.data.frame(pso_pred_FE_new$summary)
FULL_inter_prob <- rbind(FULL_inter_prob, FULL_inter_prob_new)
FULL_inter_prob <- FULL_inter_prob %>% mutate(Network = "FULL")
FULL_inter_prob <- FULL_inter_prob %>% mutate(Adjusted = "MLNMR")

# NMA
FULL_inter_probit_nma <- as.data.frame(pso_releff_FE_UA$summary)
FULL_inter_probit_nma <- duplicate_and_add_study(FULL_inter_probit_nma, study_names_int)
FULL_inter_probit_nma <- FULL_inter_probit_nma %>% mutate(Network = "FULL")
FULL_inter_probit_nma <- FULL_inter_probit_nma %>% mutate(Adjusted = "NMA")
FULL_inter_probit_nma <- add_all_categories(FULL_inter_probit_nma)

FULL_inter_prob_nma <- as.data.frame(pso_pred_FE_UA$summary)
FULL_inter_prob_nma_chir <- as.data.frame(pso_pred_FE_new_UA_chir$summary)
FULL_inter_prob_nma_chir <- FULL_inter_prob_nma_chir %>% mutate(.study = "Chiricozzi 2019")
FULL_inter_prob_nma_PROS <- as.data.frame(pso_pred_FE_new_UA_PROS$summary)
FULL_inter_prob_nma_PROS <- FULL_inter_prob_nma_PROS %>% mutate(.study = "PROSPECT")
FULL_inter_prob_nma <- rbind(FULL_inter_prob_nma, FULL_inter_prob_nma_chir, FULL_inter_prob_nma_PROS)
FULL_inter_prob_nma <- FULL_inter_prob_nma %>% mutate(Network = "FULL")
FULL_inter_prob_nma <- FULL_inter_prob_nma %>% mutate(Adjusted = "NMA")

### FIXED ----
FIXED_inter_probit <- as.data.frame(pso_releff_FE_FIXED$summary)
FIXED_inter_probit_new <- as.data.frame(pso_releff_FE_new_FIXED$summary)
FIXED_inter_probit <- rbind(FIXED_inter_probit,FIXED_inter_probit_new)
FIXED_inter_probit <- FIXED_inter_probit %>% mutate(Network = "FIXED - High ESS")
FIXED_inter_probit <- FIXED_inter_probit %>% mutate(Adjusted = "MLNMR")

FIXED_inter_prob <- as.data.frame(pso_pred_FE_FIXED$summary)
FIXED_inter_prob_new <- as.data.frame(pso_pred_FE_new_FIXED$summary)
FIXED_inter_prob <- rbind(FIXED_inter_prob, FIXED_inter_prob_new)
FIXED_inter_prob <- FIXED_inter_prob %>% mutate(Network = "FIXED - High ESS")
FIXED_inter_prob <- FIXED_inter_prob %>% mutate(Adjusted = "MLNMR")

FIXED_inter_probit_bad <- as.data.frame(pso_releff_FE_FIXED_bad$summary)
FIXED_inter_probit_new_bad <- as.data.frame(pso_releff_FE_new_FIXED_bad$summary)
FIXED_inter_probit_bad <- rbind(FIXED_inter_probit_bad,FIXED_inter_probit_new_bad)
FIXED_inter_probit_bad <- FIXED_inter_probit_bad %>% mutate(Network = "FIXED - Low ESS")
FIXED_inter_probit_bad <- FIXED_inter_probit_bad %>% mutate(Adjusted = "MLNMR")

FIXED_inter_prob_bad <- as.data.frame(pso_pred_FE_FIXED_bad$summary)
FIXED_inter_prob_new_bad <- as.data.frame(pso_pred_FE_new_FIXED_bad$summary)
FIXED_inter_prob_bad <- rbind(FIXED_inter_prob_bad, FIXED_inter_prob_new_bad)
FIXED_inter_prob_bad <- FIXED_inter_prob_bad %>% mutate(Network = "FIXED - Low ESS")
FIXED_inter_prob_bad <- FIXED_inter_prob_bad %>% mutate(Adjusted = "MLNMR")

# NMA
FIXED_inter_probit_nma <- as.data.frame(pso_releff_FE_FIXED_UA$summary)
FIXED_inter_probit_nma <- duplicate_and_add_study(FIXED_inter_probit_nma, study_names_int)
FIXED_inter_probit_nma <- FIXED_inter_probit_nma %>% mutate(Network = "FIXED - High ESS")
FIXED_inter_probit_nma <- FIXED_inter_probit_nma %>% mutate(Adjusted = "NMA")
FIXED_inter_probit_nma <- add_all_categories(FIXED_inter_probit_nma)

FIXED_inter_probit_nma_bad <- as.data.frame(pso_releff_FE_FIXED_bad_UA$summary)
FIXED_inter_probit_nma_bad <- duplicate_and_add_study(FIXED_inter_probit_nma_bad, study_names_int)
FIXED_inter_probit_nma_bad <- FIXED_inter_probit_nma_bad %>% mutate(Network = "FIXED - Low ESS")
FIXED_inter_probit_nma_bad <- FIXED_inter_probit_nma_bad %>% mutate(Adjusted = "NMA")
FIXED_inter_probit_nma_bad <- add_all_categories(FIXED_inter_probit_nma_bad)

FIXED_inter_prob_nma <- as.data.frame(pso_pred_FE_FIXED_UA$summary)
FIXED_inter_prob_nma_chir <- as.data.frame(pso_pred_FE_new_FIXED_UA_chir$summary)
FIXED_inter_prob_nma_chir <- FIXED_inter_prob_nma_chir %>% mutate(.study = "Chiricozzi 2019")
FIXED_inter_prob_nma_PROS <- as.data.frame(pso_pred_FE_new_FIXED_UA_PROS$summary)
FIXED_inter_prob_nma_PROS <- FIXED_inter_prob_nma_PROS %>% mutate(.study = "PROSPECT")
FIXED_inter_prob_nma <- rbind(FIXED_inter_prob_nma, FIXED_inter_prob_nma_chir, FIXED_inter_prob_nma_PROS)
FIXED_inter_prob_nma <- FIXED_inter_prob_nma %>% mutate(Network = "FIXED - High ESS")
FIXED_inter_prob_nma <- FIXED_inter_prob_nma %>% mutate(Adjusted = "NMA")

FIXED_inter_prob_nma_bad <- as.data.frame(pso_pred_FE_FIXED_bad_UA$summary)
FIXED_inter_prob_nma_chir_bad <- as.data.frame(pso_pred_FE_new_FIXED_bad_UA_chir$summary)
FIXED_inter_prob_nma_chir_bad <- FIXED_inter_prob_nma_chir_bad %>% mutate(.study = "Chiricozzi 2019")
FIXED_inter_prob_nma_PROS_bad <- as.data.frame(pso_pred_FE_new_FIXED_bad_UA_PROS$summary)
FIXED_inter_prob_nma_PROS_bad <- FIXED_inter_prob_nma_PROS_bad %>% mutate(.study = "PROSPECT")
FIXED_inter_prob_nma_bad <- rbind(FIXED_inter_prob_nma_bad, FIXED_inter_prob_nma_chir_bad, FIXED_inter_prob_nma_PROS_bad)
FIXED_inter_prob_nma_bad <- FIXED_inter_prob_nma_bad %>% mutate(Network = "FIXED - Low ESS")
FIXED_inter_prob_nma_bad <- FIXED_inter_prob_nma_bad %>% mutate(Adjusted = "NMA")

### RANDOM -----
### 1P ----
RANDOM_1P_inter_probit <- as.data.frame(pso_releff_FE_RANDOM_1P$summary)
RANDOM_1P_inter_probit_new <- as.data.frame(pso_releff_FE_new_RANDOM_1P$summary)
RANDOM_1P_inter_probit <- rbind(RANDOM_1P_inter_probit, RANDOM_1P_inter_probit_new)
RANDOM_1P_inter_probit <- RANDOM_1P_inter_probit %>% mutate(Network = "RANDOM_1P")
RANDOM_1P_inter_probit <- RANDOM_1P_inter_probit %>% mutate(Adjusted = "MLNMR")

RANDOM_1P_inter_prob <- as.data.frame(pso_pred_FE_RANDOM_1P$summary)
RANDOM_1P_inter_prob_new <- as.data.frame(pso_pred_FE_new_RANDOM_1P$summary)
RANDOM_1P_inter_prob <- rbind(RANDOM_1P_inter_prob, RANDOM_1P_inter_prob_new)
RANDOM_1P_inter_prob <- RANDOM_1P_inter_prob %>% mutate(Network = "RANDOM_1P")
RANDOM_1P_inter_prob <- RANDOM_1P_inter_prob %>% mutate(Adjusted = "MLNMR")

# NMA
RANDOM_1P_inter_probit_nma <- as.data.frame(pso_releff_FE_RANDOM_UA_1P$summary)
RANDOM_1P_inter_probit_nma <- duplicate_and_add_study(RANDOM_1P_inter_probit_nma, study_names_int)
RANDOM_1P_inter_probit_nma <- RANDOM_1P_inter_probit_nma %>% mutate(Network = "RANDOM_1P")
RANDOM_1P_inter_probit_nma <- RANDOM_1P_inter_probit_nma %>% mutate(Adjusted = "NMA")
RANDOM_1P_inter_probit_nma <- add_all_categories(RANDOM_1P_inter_probit_nma)

RANDOM_1P_inter_prob_nma <- as.data.frame(pso_pred_FE_RANDOM_UA_1P$summary)
RANDOM_1P_inter_prob_nma_chir <- as.data.frame(pso_pred_FE_new_RANDOM_UA_chir_1P$summary)
RANDOM_1P_inter_prob_nma_chir <- RANDOM_1P_inter_prob_nma_chir %>% mutate(.study = "Chiricozzi 2019")
RANDOM_1P_inter_prob_nma_PROS <- as.data.frame(pso_pred_FE_new_RANDOM_UA_PROS_1P$summary)
RANDOM_1P_inter_prob_nma_PROS <- RANDOM_1P_inter_prob_nma_PROS %>% mutate(.study = "PROSPECT")
RANDOM_1P_inter_prob_nma <- rbind(RANDOM_1P_inter_prob_nma, RANDOM_1P_inter_prob_nma_chir, RANDOM_1P_inter_prob_nma_PROS)
RANDOM_1P_inter_prob_nma <- RANDOM_1P_inter_prob_nma %>% mutate(Network = "RANDOM_1P")
RANDOM_1P_inter_prob_nma <- RANDOM_1P_inter_prob_nma %>% mutate(Adjusted = "NMA")

### 3P ----
RANDOM_3P_inter_probit <- as.data.frame(pso_releff_FE_RANDOM_3P$summary)
RANDOM_3P_inter_probit_new <- as.data.frame(pso_releff_FE_new_RANDOM_3P$summary)
RANDOM_3P_inter_probit <- rbind(RANDOM_3P_inter_probit, RANDOM_3P_inter_probit_new)
RANDOM_3P_inter_probit <- RANDOM_3P_inter_probit %>% mutate(Network = "RANDOM_3P")
RANDOM_3P_inter_probit <- RANDOM_3P_inter_probit %>% mutate(Adjusted = "MLNMR")
RANDOM_3P_inter_probit <- RANDOM_3P_inter_probit
RANDOM_3P_inter_prob <- as.data.frame(pso_pred_FE_RANDOM_3P$summary)
RANDOM_3P_inter_prob_new <- as.data.frame(pso_pred_FE_new_RANDOM_3P$summary)
RANDOM_3P_inter_prob <- rbind(RANDOM_3P_inter_prob, RANDOM_3P_inter_prob_new)
RANDOM_3P_inter_prob <- RANDOM_3P_inter_prob %>% mutate(Network = "RANDOM_3P")
RANDOM_3P_inter_prob <- RANDOM_3P_inter_prob %>% mutate(Adjusted = "MLNMR")

# NMA
RANDOM_3P_inter_probit_nma <- as.data.frame(pso_releff_FE_RANDOM_UA_3P$summary)
RANDOM_3P_inter_probit_nma <- duplicate_and_add_study(RANDOM_3P_inter_probit_nma, study_names_int)
RANDOM_3P_inter_probit_nma <- RANDOM_3P_inter_probit_nma %>% mutate(Network = "RANDOM_3P")
RANDOM_3P_inter_probit_nma <- RANDOM_3P_inter_probit_nma %>% mutate(Adjusted = "NMA")
RANDOM_3P_inter_probit_nma <- add_all_categories(RANDOM_3P_inter_probit_nma)

RANDOM_3P_inter_prob_nma <- as.data.frame(pso_pred_FE_RANDOM_UA_3P$summary)
RANDOM_3P_inter_prob_nma_chir <- as.data.frame(pso_pred_FE_new_RANDOM_UA_chir_3P$summary)
RANDOM_3P_inter_prob_nma_chir <- RANDOM_3P_inter_prob_nma_chir %>% mutate(.study = "Chiricozzi 2019")
RANDOM_3P_inter_prob_nma_PROS <- as.data.frame(pso_pred_FE_new_RANDOM_UA_PROS_3P$summary)
RANDOM_3P_inter_prob_nma_PROS <- RANDOM_3P_inter_prob_nma_PROS %>% mutate(.study = "PROSPECT")
RANDOM_3P_inter_prob_nma <- rbind(RANDOM_3P_inter_prob_nma, RANDOM_3P_inter_prob_nma_chir, RANDOM_3P_inter_prob_nma_PROS)
RANDOM_3P_inter_prob_nma <- RANDOM_3P_inter_prob_nma %>% mutate(Network = "RANDOM_3P")
RANDOM_3P_inter_prob_nma <- RANDOM_3P_inter_prob_nma %>% mutate(Adjusted = "NMA")


# External control PROS
External_PROS_inter_probit <- as.data.frame(pso_releff_FE_external_PROS$summary)
External_PROS_inter_probit_new <- as.data.frame(pso_releff_FE_new_external_PROS)
External_PROS_inter_probit <- rbind(External_PROS_inter_probit, External_PROS_inter_probit_new)
External_PROS_inter_probit <- External_PROS_inter_probit %>% mutate(Network = "EXTERNAL_Pros")
External_PROS_inter_probit <- External_PROS_inter_probit %>% mutate(Adjusted = "MLNMR")

External_PROS_inter_prob <- as.data.frame(pso_pred_FE_external_PROS$summary)
External_PROS_inter_prob_new <- as.data.frame(pso_pred_FE_new_external_PROS$summary)
External_PROS_inter_prob <- rbind(External_PROS_inter_prob, External_PROS_inter_prob_new)
External_PROS_inter_prob <- External_PROS_inter_prob %>% mutate(Network = "EXTERNAL_Pros")
External_PROS_inter_prob <- External_PROS_inter_prob %>% mutate(Adjusted = "MLNMR")

# NMA
External_PROS_inter_probit_nma <- as.data.frame(pso_releff_FE_external_PROS_UA$summary)
External_PROS_inter_probit_nma <- duplicate_and_add_study(External_PROS_inter_probit_nma, study_names_int)
External_PROS_inter_probit_nma <- External_PROS_inter_probit_nma %>% mutate(Network = "EXTERNAL_Pros")
External_PROS_inter_probit_nma <- External_PROS_inter_probit_nma %>% mutate(Adjusted = "NMA")
External_PROS_inter_probit_nma <- add_all_categories(External_PROS_inter_probit_nma)

External_PROS_inter_prob_nma <- as.data.frame(pso_pred_FE_external_PROS_UA$summary)
External_PROS_inter_prob_nma_chir <- as.data.frame(pso_pred_FE_new_external_PROS_UA_chir$summary)
External_PROS_inter_prob_nma_chir <- External_PROS_inter_prob_nma_chir %>% mutate(.study = "Chiricozzi 2019")
External_PROS_inter_prob_nma_PROS <- as.data.frame(pso_pred_FE_new_external_PROS_UA_PROS$summary)
External_PROS_inter_prob_nma_PROS <- External_PROS_inter_prob_nma_PROS %>% mutate(.study = "PROSPECT")
External_PROS_inter_prob_nma <- rbind(External_PROS_inter_prob_nma, External_PROS_inter_prob_nma_chir, External_PROS_inter_prob_nma_PROS)
External_PROS_inter_prob_nma <- External_PROS_inter_prob_nma %>% mutate(Network = "EXTERNAL_Pros")
External_PROS_inter_prob_nma <- External_PROS_inter_prob_nma %>% mutate(Adjusted = "NMA")

# External control chir
External_chir_inter_probit <- as.data.frame(pso_releff_FE_external_chir$summary)
External_chir_inter_probit_new <- as.data.frame(pso_releff_FE_new_external_chir$summary)
External_chir_inter_probit <- rbind(External_chir_inter_probit, External_chir_inter_probit_new)
External_chir_inter_probit <- External_chir_inter_probit %>% mutate(Network = "EXTERNAL_Chir")
External_chir_inter_probit <- External_chir_inter_probit %>% mutate(Adjusted = "MLNMR")

External_chir_inter_prob <- as.data.frame(pso_pred_FE_external_chir$summary)
External_chir_inter_prob_new <- as.data.frame(pso_pred_FE_new_external_chir$summary)
External_chir_inter_prob <- rbind(External_chir_inter_prob, External_chir_inter_prob_new)
External_chir_inter_prob <- External_chir_inter_prob %>% mutate(Network = "EXTERNAL_Chir")
External_chir_inter_prob <- External_chir_inter_prob %>% mutate(Adjusted = "MLNMR")

# NMA
External_chir_inter_probit_nma <- as.data.frame(pso_releff_FE_external_chir_UA$summary)
External_chir_inter_probit_nma <- duplicate_and_add_study(External_chir_inter_probit_nma, study_names_int)
External_chir_inter_probit_nma <- External_chir_inter_probit_nma %>% mutate(Network = "EXTERNAL_Chir")
External_chir_inter_probit_nma <- External_chir_inter_probit_nma %>% mutate(Adjusted = "NMA")
External_chir_inter_probit_nma <- add_all_categories(External_chir_inter_probit_nma)

External_chir_inter_prob_nma <- as.data.frame(pso_pred_FE_external_chir_UA$summary)
External_chir_inter_prob_nma_chir <- as.data.frame(pso_pred_FE_new_external_chir_UA_chir$summary)
External_chir_inter_prob_nma_chir <- External_chir_inter_prob_nma_chir %>% mutate(.study = "Chiricozzi 2019")
External_chir_inter_prob_nma_PROS <- as.data.frame(pso_pred_FE_new_external_chir_UA_PROS$summary)
External_chir_inter_prob_nma_PROS <- External_chir_inter_prob_nma_PROS %>% mutate(.study = "PROSPECT")
External_chir_inter_prob_nma <- rbind(External_chir_inter_prob_nma, External_chir_inter_prob_nma_chir, External_chir_inter_prob_nma_PROS)
External_chir_inter_prob_nma <- External_chir_inter_prob_nma %>% mutate(Network = "EXTERNAL_Chir")
External_chir_inter_prob_nma <- External_chir_inter_prob_nma %>% mutate(Adjusted = "NMA")

ALL_prob <- rbind(FULL_inter_prob,
                  FULL_inter_prob_nma,
                  FIXED_inter_prob,
                  FIXED_inter_prob_nma,
                  FIXED_inter_prob_bad,
                  FIXED_inter_prob_nma_bad,
                  RANDOM_1P_inter_prob,
                  RANDOM_1P_inter_prob_nma,
                  RANDOM_3P_inter_prob,
                  RANDOM_3P_inter_prob_nma,
                  External_PROS_inter_prob,
                  External_PROS_inter_prob_nma,
                  External_chir_inter_prob,
                  External_chir_inter_prob_nma)

ALL_prob <- ALL_prob %>% filter(.category %in% pasi_categories)

ALL_probit <- rbind(FULL_inter_probit,
                    FULL_inter_probit_nma,
                    FIXED_inter_probit,
                    FIXED_inter_probit_nma,
                    FIXED_inter_probit_bad,
                    FIXED_inter_probit_nma_bad,
                    RANDOM_1P_inter_probit,
                    RANDOM_1P_inter_probit_nma,
                    RANDOM_3P_inter_probit,
                    RANDOM_3P_inter_probit_nma,
                    External_PROS_inter_probit,
                    External_PROS_inter_probit_nma,
                    External_chir_inter_probit,
                    External_chir_inter_probit_nma)

ALL_prob$.category <- factor(ALL_prob$.category, levels = pasi_categories)
ALL_probit$.category <- factor(ALL_probit$.category, levels = pasi_categories)


p_trt_label <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Treatment", angle = -90, size = 11 / .pt) +
  theme_void()

recode_network <- function(data) {
  data %>%
    mutate(
      Network = as.character(Network),
      Network = case_when(
        Network == "FULL"             ~ "**Full**",
        Network == "FIXED - High ESS" ~ "Fixed (UNC-3/ERA)",
        Network == "FIXED - Low ESS"  ~ "Fixed (JUNC/IXOR)",
        Network == "RANDOM_1P"        ~ "Random (1P)",
        Network == "RANDOM_3P"        ~ "Random (3P)",
        Network == "EXTERNAL_Pros"    ~ "External (PROSPECT)",
        Network == "EXTERNAL_Chir"    ~ "External (Chiricozzi)",
        TRUE                          ~ Network
      )
    )
}

# Probability Plots ----
# One plot per PASI threshold: PASI75 is the main-text figure, PASI90/PASI100
# are the appendix figures.

make_prob_plot <- function(category, title = paste("Target population -", category)) {
  FILTERED <- ALL_prob %>%
    filter(.category == category) %>%
    filter(!(.study %in% c("UNCOVER-3 & ERASURE", "JUNCTURE & IXORA-S", "PROSPECT & IXORA-S", "Chiricozzi 2019 & IXORA-S", "Chiricozzi 2019", "UNCOVER-1", "UNCOVER-2", "UNCOVER-3", "JUNCTURE", "FIXTURE & ERASURE", "CLEAR", "ERASURE", "IXORA-S"))) %>%
    recode_network()

  p_prob <- FILTERED %>%
    rename(Treatment = .trt, `Target population` = .study) %>%
    ggplot(aes(x=Network, y=mean, ymin=`2.5%`, ymax=`97.5%`, color = Adjusted)) +
    geom_point(position = position_dodge(width = .5), size = 0.1) +
    geom_pointrange(position = position_dodge(width = .5), size = 0.1) +
    coord_flip() +
    xlab("Bridge") + ylab("Probability") +
    scale_y_continuous(breaks = c(0.25, 0.50, 0.75, 1.00)) +
    labs(color = "Model") +
    facet_grid(Treatment ~ `Target population`, scales = "free_y",
               labeller = label_value) +
    ggtitle(title) +
    theme_bw() +
    theme(panel.spacing = unit(0.5, "lines"),
          plot.title = element_text(size = 11, hjust = 0.5),
          axis.title = element_text(size = 11),

          strip.text.y = element_text(angle = 0, hjust = 0, vjust = 1),
          axis.text.y = element_markdown(hjust = 1, margin = margin(r = 3, unit = "mm")),
          legend.position = c(0.97, 0.97),
          legend.justification = c(1, 1),
          legend.background = element_rect(fill = "white", colour = "grey80"),
          legend.margin = margin(4, 4, 4, 4))

  (p_prob + p_trt_label) + plot_layout(widths = c(25, 1))
}

# Main text
FILTERED_prob_plot_75 <- make_prob_plot("PASI75")
FILTERED_prob_plot_75

# Appendix
FILTERED_prob_plot_90 <- make_prob_plot("PASI90")
FILTERED_prob_plot_90

FILTERED_prob_plot_100 <- make_prob_plot("PASI100")
FILTERED_prob_plot_100

# Probit Plots ----
make_probit_plot <- function(category, title = paste("Target population -", category)) {
  FILTERED_probit <- ALL_probit %>%
    filter(.category == category) %>%
    filter(!(.study %in% c("UNCOVER-3 & ERASURE", "PsoBest", "JUNCTURE & IXORA-S", "PROSPECT & IXORA-S", "Chiricozzi 2019 & IXORA-S", "Chiricozzi 2019", "UNCOVER-1", "UNCOVER-2", "UNCOVER-3", "JUNCTURE", "FIXTURE & ERASURE", "CLEAR", "ERASURE", "IXORA-S"))) %>%
    recode_network()

  p_probit <- FILTERED_probit %>%
    rename(Treatment = .trtb, `Target population` = .study) %>%
    ggplot(aes(x=Network, y=mean, ymin=`2.5%`, ymax=`97.5%`, color = Adjusted)) +
    geom_point(position = position_dodge(width = .5), size = 0.1) +
    geom_pointrange(position = position_dodge(width = .5), size = 0.1) +
    coord_flip() +
    xlab("Bridge") + ylab("Marginal probit difference") +
    labs(color = "Model") +
    facet_grid(Treatment ~ `Target population`, scales = "free_y",
               labeller = label_value) +
    ggtitle(title) +
    theme_bw() +
    theme(panel.spacing = unit(0.5, "lines"),
          plot.title = element_text(size = 11, hjust = 0.5),
          axis.title = element_text(size = 11),

          strip.text.y = element_text(angle = 0, hjust = 0, vjust = 1),
          axis.text.y = element_markdown(hjust = 1, margin = margin(r = 3, unit = "mm")),
          legend.position = c(0.97, 0.97),
          legend.justification = c(1, 1),
          legend.background = element_rect(fill = "white", colour = "grey80"),
          legend.margin = margin(4, 4, 4, 4))

  (p_probit + p_trt_label) + plot_layout(widths = c(25, 1))
}

# Main text
FILTERED_probit_plot_75 <- make_probit_plot("PASI75")
FILTERED_probit_plot_75

# Appendix
FILTERED_probit_plot_90 <- make_probit_plot("PASI90")
FILTERED_probit_plot_90

FILTERED_probit_plot_100 <- make_probit_plot("PASI100")
FILTERED_probit_plot_100

# Save plots ----
plots_dir <- "Plots"
dir.create(plots_dir, showWarnings = FALSE)

ggsave(file.path(plots_dir, "Probability_plot_PASI75.pdf"), FILTERED_prob_plot_75,
       width = 10, height = 12, units = "in", dpi = 300)
ggsave(file.path(plots_dir, "Probability_plot_PASI90.pdf"), FILTERED_prob_plot_90,
       width = 10, height = 12, units = "in", dpi = 300)
ggsave(file.path(plots_dir, "Probability_plot_PASI100.pdf"), FILTERED_prob_plot_100,
       width = 10, height = 12, units = "in", dpi = 300)

ggsave(file.path(plots_dir, "Probit_plot_PASI75.pdf"), FILTERED_probit_plot_75,
       width = 10, height = 12, units = "in", dpi = 300)
ggsave(file.path(plots_dir, "Probit_plot_PASI90.pdf"), FILTERED_probit_plot_90,
       width = 10, height = 12, units = "in", dpi = 300)
ggsave(file.path(plots_dir, "Probit_plot_PASI100.pdf"), FILTERED_probit_plot_100,
       width = 10, height = 12, units = "in", dpi = 300)
