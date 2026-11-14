# =====================================================================
# Marine reserve DCE - robustness checks for JEMECP-D-26-00394
# RPL-EC in utility space, CORRELATED parameters, gender-weighted (Apollo)
# Same specification as the NLOGIT RPL-EC in the paper (Tables 3, 4, 8):
#   U(one), U(two) = MR*mr + Job*job + Carbon*carbon + Cost*cost + EC
#   U(neither)     = b0
#   b0, MR, Job, Carbon ~ normal; Cost ~ triangular; all five correlated
#   (Cholesky, 15 elements); EC shared by one and two; 21 parameters.
#
# Runs three samples in sequence and builds the comparison table:
#   baseline   : exclall3 == 0                  (371)
#   panel      : baseline, online panel only    (312)
#   protest_in : drops only Excl_sca == 1       (435)
# =====================================================================
setwd("C:/Users/mp1149/Documents/12-2020/Marine Project/JEM submission/reviewer's comments/DATA_Caplan_20260116_ALL/DATA_Caplan_20260116_ALL")
rm(list = ls())
library(apollo)

N_DRAWS  <- 1000                     # as in NLOGIT (use 500 to go faster)
SAMPLES  <- c("baseline", "panel", "protest_in")
DATAFILE <- "MR_dataset.csv"

# ---------------------------------------------------------------------
# Paper values: RPL-EC (utility space, correlated), Table 4
# ---------------------------------------------------------------------
paper <- c(wtp_sq = -168.329, wtp_mr = 3.213, wtp_job = 0.670, wtp_carb = 1.735)
paper <- c(paper,
  optimistic  = unname(50 * paper["wtp_mr"] + 200 * paper["wtp_job"] + 100 * paper["wtp_carb"]),
  moderate    = unname(25 * paper["wtp_mr"] + 100 * paper["wtp_job"] +  50 * paper["wtp_carb"]),
  pessimistic = unname(-50 * paper["wtp_mr"] - 200 * paper["wtp_job"]))
paper_se <- c(wtp_sq = 37.332, wtp_mr = 0.396, wtp_job = 0.103, wtp_carb = 0.367,
              optimistic = NA, moderate = NA, pessimistic = NA)   # fill scenario SEs from Table 8

raw0 <- read.csv(DATAFILE)
raw0 <- raw0[order(raw0$pid, raw0$gid, raw0$alt), ]
resp <- raw0[!duplicated(raw0$pid), ]

results <- list()

for (SAMPLE in SAMPLES) {

  # -------------------------------------------------------------------
  # 1. Data
  # -------------------------------------------------------------------
  keep <- switch(SAMPLE,
    baseline   = resp$pid[resp$exclall3 == 0],
    panel      = resp$pid[resp$exclall3 == 0 & resp$Wave != 4],
    protest_in = resp$pid[resp$Excl_sca == 0])
  raw <- raw0[raw0$pid %in% keep, ]
  w <- reshape(raw[, c("pid", "gid", "alt", "Y", "mr", "job", "carbon", "cost")],
               idvar = c("pid", "gid"), timevar = "alt", direction = "wide", sep = "_")
  w$choice <- 1 * w$Y_1 + 2 * w$Y_2 + 3 * w$Y_3
  wt <- resp[match(w$pid, resp$pid), "wtgend"]
  w$wtgend <- wt / mean(resp$wtgend[resp$pid %in% keep])
  database <- w[order(w$pid, w$gid), ]
  cat("\n=====", SAMPLE, "| respondents:", length(unique(database$pid)), "=====\n")

  # -------------------------------------------------------------------
  # 2. Model
  # -------------------------------------------------------------------
  apollo_initialise()
  apollo_control <- list(
    modelName  = paste0("RPLEC_corr_", SAMPLE),
    modelDescr = "RPL-EC utility space, correlated, gender weights",
    indivID    = "pid", mixing = TRUE, weights = "wtgend",
    nCores     = 7, analyticGrad = FALSE, noModification = TRUE,
    outputDirectory = "output")

  # start values: paper means; Cholesky diagonal from paper SDs, off-diagonals 0
  apollo_beta <- c(mu_sq = -1.003, mu_mr = 0.019, mu_job = 0.004, mu_carb = 0.010, mu_cost = -0.006,
                    l11 = 2.101,
                    l21 = 0, l22 = 0.019,
                    l31 = 0, l32 = 0, l33 = 0.003,
                    l41 = 0, l42 = 0, l43 = 0, l44 = 0.010,
                    l51 = 0, l52 = 0, l53 = 0, l54 = 0, l55 = 0.012,
                    sig_ec = 0.462)
  apollo_fixed <- c()

  apollo_draws <- list(
    interDrawsType = "halton", interNDraws = N_DRAWS,
    interUnifDraws = c("u_c1", "u_c2"),
    interNormDraws = c("v1", "v2", "v3", "v4", "d_ec"))

  apollo_randCoeff <- function(apollo_beta, apollo_inputs) {
    rc <- list()
    v5 <- u_c1 + u_c2 - 1                                   # triangular(-1,1) for cost
    rc[["b_sq"]]   <- mu_sq   + l11 * v1
    rc[["b_mr"]]   <- mu_mr   + l21 * v1 + l22 * v2
    rc[["b_job"]]  <- mu_job  + l31 * v1 + l32 * v2 + l33 * v3
    rc[["b_carb"]] <- mu_carb + l41 * v1 + l42 * v2 + l43 * v3 + l44 * v4
    rc[["b_cost"]] <- mu_cost + l51 * v1 + l52 * v2 + l53 * v3 + l54 * v4 + l55 * v5
    rc[["ec"]]     <- sig_ec * d_ec
    return(rc)
  }

  apollo_inputs <- apollo_validateInputs()

  apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    P <- list(); V <- list()
    V[["alt1"]] <- b_mr * mr_1 + b_job * job_1 + b_carb * carbon_1 + b_cost * cost_1 + ec
    V[["alt2"]] <- b_mr * mr_2 + b_job * job_2 + b_carb * carbon_2 + b_cost * cost_2 + ec
    V[["alt3"]] <- b_sq
    mnl_settings <- list(alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
                         avail = list(alt1 = 1, alt2 = 1, alt3 = 1),
                         choiceVar = choice, utilities = V)
    P[["model"]] <- apollo_mnl(mnl_settings, functionality)
    P <- apollo_panelProd(P, apollo_inputs, functionality)
    P <- apollo_avgInterDraws(P, apollo_inputs, functionality)
    P <- apollo_weighting(P, apollo_inputs, functionality)
    P <- apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
  }

  model <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)
  apollo_modelOutput(model)

  # -------------------------------------------------------------------
  # 3. WTP and SDs (delta method, robust SEs)
  # -------------------------------------------------------------------
  dm <- apollo_deltaMethod(model, list(expression = c(
    wtp_sq      = "-mu_sq/mu_cost",
    wtp_mr      = "-mu_mr/mu_cost",
    wtp_job     = "-mu_job/mu_cost",
    wtp_carb    = "-mu_carb/mu_cost",
    optimistic  = "-(50*mu_mr + 200*mu_job + 100*mu_carb)/mu_cost",
    moderate    = "-(25*mu_mr + 100*mu_job +  50*mu_carb)/mu_cost",
    pessimistic = "-(-50*mu_mr - 200*mu_job)/mu_cost",
    sd_mr       = "sqrt(l21^2 + l22^2)",
    sd_job      = "sqrt(l31^2 + l32^2 + l33^2)",
    sd_carb     = "sqrt(l41^2 + l42^2 + l43^2 + l44^2)")))
  results[[SAMPLE]] <- list(dm = dm, LL = model$maximum, N = length(unique(database$pid)))
}

# ---------------------------------------------------------------------
# 4. Comparison table: paper | Apollo baseline | panel | protest_in
#    % changes relative to the Apollo baseline (same software and model)
# ---------------------------------------------------------------------
fmt <- function(v, s) ifelse(is.na(v), "", ifelse(is.na(s), sprintf("%.2f", v), sprintf("%.2f (%.2f)", v, s)))
rows <- c("wtp_sq", "wtp_mr", "wtp_job", "wtp_carb", "optimistic", "moderate", "pessimistic",
          "sd_mr", "sd_job", "sd_carb")
get <- function(s, col) setNames(results[[s]]$dm[[col]], results[[s]]$dm$Expression)[rows]
b  <- get("baseline", "Value");   bs <- get("baseline", "s.e.")
pn <- get("panel", "Value");      ps <- get("panel", "s.e.")
pr <- get("protest_in", "Value"); rs <- get("protest_in", "s.e.")
pct <- function(x) sprintf("%+.0f%%", 100 * (x / b - 1))

tab <- data.frame(
  Row               = rows,
  Paper_NLOGIT      = fmt(paper[rows], paper_se[rows]),
  Apollo_baseline   = fmt(b, bs),
  Panel_only        = paste(fmt(pn, ps), pct(pn)),
  Protesters_in     = paste(fmt(pr, rs), pct(pr)))
tab <- rbind(tab, data.frame(Row = c("LL", "N"),
  Paper_NLOGIT    = c("-1299.93", "371"),
  Apollo_baseline = c(sprintf("%.2f", results$baseline$LL),   results$baseline$N),
  Panel_only      = c(sprintf("%.2f", results$panel$LL),      results$panel$N),
  Protesters_in   = c(sprintf("%.2f", results$protest_in$LL), results$protest_in$N)))
print(tab, row.names = FALSE)
write.csv(tab, "output/robustness_table_corr.csv", row.names = FALSE)
