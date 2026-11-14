# =====================================================================
# Marine reserve DCE - non-linear RPL-EC in WTP space + Wald tests (Apollo)
# For R2 #4 (gain-loss asymmetry) on the preferred WTP-space specification
#
# Same structure as the NLOGIT non-linear RPL-EC 2 (WTP space), Tables 5-6:
#   U(one), U(two) = a * ( -w_mr1*MR(-50) - w_mr3*MR(+50) - w_jb1*Jobs(-200)
#                          - w_jb3*Jobs(+200) - w_c2*C(+25) - w_c3*C(+50)
#                          - w_c4*C(+100) + Cost ) + EC
#   U(neither)     = a * ( -b )
#   b and the 7 level WTPs ~ normal (uncorrelated); a ~ triangular (free spread);
#   EC shared by options 1 and 2; gender weights; panel over 4 tasks.
#   Base levels: MR 0%, Jobs 0, Carbon 0% (status quo).
# In WTP space the level coefficients ARE the WTPs (US$ per household per year).
# =====================================================================
setwd("C:/Users/mp1149/Documents/12-2020/Marine Project/JEM submission/reviewer's comments/DATA_Caplan_20260116_ALL/DATA_Caplan_20260116_ALL")
rm(list = ls())
library(apollo)

N_DRAWS  <- 1000            # as in NLOGIT (500 to go faster)
DATAFILE <- "MR_dataset.csv"

# ---------------------------------------------------------------------
# 1. Data: paper sample, long -> wide, level dummies
# ---------------------------------------------------------------------
raw  <- read.csv(DATAFILE)
raw  <- raw[order(raw$pid, raw$gid, raw$alt), ]
resp <- raw[!duplicated(raw$pid), ]
keep <- resp$pid[resp$exclall3 == 0]
raw  <- raw[raw$pid %in% keep, ]

raw$mr1 <- 1 * (raw$mr  == -50);  raw$mr3 <- 1 * (raw$mr  ==  50)
raw$jb1 <- 1 * (raw$job == -200); raw$jb3 <- 1 * (raw$job == 200)
raw$c2  <- 1 * (raw$carbon == 25); raw$c3 <- 1 * (raw$carbon == 50); raw$c4 <- 1 * (raw$carbon == 100)

w <- reshape(raw[, c("pid", "gid", "alt", "Y", "mr1", "mr3", "jb1", "jb3", "c2", "c3", "c4", "cost")],
             idvar = c("pid", "gid"), timevar = "alt", direction = "wide", sep = "_")
w$choice <- 1 * w$Y_1 + 2 * w$Y_2 + 3 * w$Y_3
w$wtgend <- resp$wtgend[match(w$pid, resp$pid)] / mean(resp$wtgend[resp$pid %in% keep])
database <- w[order(w$pid, w$gid), ]
cat("Respondents:", length(unique(database$pid)), "| tasks:", nrow(database), "\n")

# ---------------------------------------------------------------------
# 2. Model
# ---------------------------------------------------------------------
apollo_initialise()
apollo_control <- list(
  modelName  = "RPLEC_WTP_nonlinear",
  modelDescr = "Non-linear RPL-EC, WTP space, uncorrelated, gender weights",
  indivID    = "pid", mixing = TRUE, weights = "wtgend",
  nCores     = 7, analyticGrad = FALSE, noModification = TRUE,
  outputDirectory = "output")

# start values: NLOGIT estimates of the non-linear RPL-EC 2 (WTP space), Apr 01 2025 output
apollo_beta <- c(mu_b  = -155.030, w_mr1 = -204.357, w_mr3 = 87.119,  w_jb1 = -139.820, w_jb3 = 138.141,
                 w_c2  = 277.724,  w_c3  = 275.671,  w_c4  = 383.165, mu_a  = -0.00371,
                 s_b   = 12.731,   s_mr1 = 225.398,  s_mr3 = 204.593, s_jb1 = 121.234,  s_jb3 = 133.628,
                 s_c2  = 504.436,  s_c3  = 323.894,  s_c4  = 702.967, s_a   = 0.00218,
                 sig_ec = 1.492)
apollo_fixed <- c()

apollo_draws <- list(
  interDrawsType = "halton", interNDraws = N_DRAWS,
  interUnifDraws = c("u_a1", "u_a2"),
  interNormDraws = c("d_b", "d_mr1", "d_mr3", "d_jb1", "d_jb3", "d_c2", "d_c3", "d_c4", "d_ec"))

apollo_randCoeff <- function(apollo_beta, apollo_inputs) {
  rc  <- list()
  tri <- u_a1 + u_a2 - 1                        # triangular(-1,1)
  rc[["b"]]   <- mu_b  + s_b   * d_b
  rc[["v1"]]  <- w_mr1 + s_mr1 * d_mr1
  rc[["v2"]]  <- w_mr3 + s_mr3 * d_mr3
  rc[["v3"]]  <- w_jb1 + s_jb1 * d_jb1
  rc[["v4"]]  <- w_jb3 + s_jb3 * d_jb3
  rc[["v5"]]  <- w_c2  + s_c2  * d_c2
  rc[["v6"]]  <- w_c3  + s_c3  * d_c3
  rc[["v7"]]  <- w_c4  + s_c4  * d_c4
  rc[["a"]]   <- mu_a  + s_a   * tri            # free spread, as NLOGIT a(t,1)
  rc[["ec"]]  <- sig_ec * d_ec
  return(rc)
}

apollo_inputs <- apollo_validateInputs()

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P <- list(); V <- list()
  V[["alt1"]] <- a * (-v1 * mr1_1 - v2 * mr3_1 - v3 * jb1_1 - v4 * jb3_1 -
                       v5 * c2_1 - v6 * c3_1 - v7 * c4_1 + cost_1) + ec
  V[["alt2"]] <- a * (-v1 * mr1_2 - v2 * mr3_2 - v3 * jb1_2 - v4 * jb3_2 -
                       v5 * c2_2 - v6 * c3_2 - v7 * c4_2 + cost_2) + ec
  V[["alt3"]] <- a * (-b)
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
apollo_modelOutput(model)   # w_* means = WTP by level; compare with Table 6, WTP-space column

# ---------------------------------------------------------------------
# 3. Wald tests (delta method, robust SEs)
# ---------------------------------------------------------------------
# Gain-loss asymmetry (R2 #4): H0: |WTP(loss)| = WTP(gain)  <=>  WTP(loss) + WTP(gain) = 0.
# A significantly NEGATIVE value means losses weigh more than equal-sized gains.
cat("\n===== Gain-loss asymmetry tests (H0: value = 0) =====\n")
apollo_deltaMethod(model, list(expression = c(
  asym_mr   = "w_mr1 + w_mr3",       # WTP(MR -50%) + WTP(MR +50%)
  asym_jobs = "w_jb1 + w_jb3")))     # WTP(Jobs -200) + WTP(Jobs +200)

# Differences between carbon levels (R1 #1): H0: no difference in WTP
cat("\n===== Differences between carbon levels (H0: value = 0) =====\n")
apollo_deltaMethod(model, list(expression = c(
  c50_minus_c25  = "w_c3 - w_c2",
  c100_minus_c50 = "w_c4 - w_c3",
  c100_minus_c25 = "w_c4 - w_c2")))
