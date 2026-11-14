# =====================================================================
# Marine reserve DCE - non-linear RPL-EC (utility space) + Wald tests
# For reviewer comments R2 #4 (gain-loss asymmetry) and R1 #1 (carbon levels)
#
# Same structure as the NLOGIT non-linear RPL-EC (Tables 5-6), but with
# UNCORRELATED random parameters:
#   U(one), U(two) = b0 + MR(-50) + MR(+50) + Jobs(-200) + Jobs(+200)
#                    + Carbon(+25) + Carbon(+50) + Carbon(+100) + Cost*cost + EC
#   U(neither)     = 0
#   b0 and the 7 level dummies ~ normal; Cost ~ triangular (free spread);
#   EC shared by options 1 and 2; gender weights; panel over 4 tasks.
#   Base levels: MR 0%, Jobs 0, Carbon 0% (status quo).
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
  modelName  = "RPLEC_nonlinear_uncorr",
  modelDescr = "Non-linear RPL-EC, utility space, uncorrelated, gender weights",
  indivID    = "pid", mixing = TRUE, weights = "wtgend",
  nCores     = 7, analyticGrad = FALSE, noModification = TRUE,
  outputDirectory = "output")

# start values: about twice the NLOGIT MNL estimates (start-value model in the output)
apollo_beta <- c(mu_b0 = -0.2, mu_mr1 = -1.1, mu_mr3 = 0.6, mu_jb1 = -0.6, mu_jb3 = 0.8,
                 mu_c2 = 1.5,  mu_c3  = 1.7,  mu_c4  = 1.9, mu_cost = -0.0045,
                 s_b0  = 1,    s_mr1  = 1,    s_mr3  = 1,   s_jb1  = 1,    s_jb3 = 1,
                 s_c2  = 1,    s_c3   = 1,    s_c4   = 1,   s_cost = 0.004,
                 sig_ec = 1)
apollo_fixed <- c()

apollo_draws <- list(
  interDrawsType = "halton", interNDraws = N_DRAWS,
  interUnifDraws = c("u_c1", "u_c2"),
  interNormDraws = c("d_b0", "d_mr1", "d_mr3", "d_jb1", "d_jb3", "d_c2", "d_c3", "d_c4", "d_ec"))

apollo_randCoeff <- function(apollo_beta, apollo_inputs) {
  rc  <- list()
  tri <- u_c1 + u_c2 - 1                       # triangular(-1,1)
  rc[["b0"]]   <- mu_b0  + s_b0  * d_b0
  rc[["bmr1"]] <- mu_mr1 + s_mr1 * d_mr1
  rc[["bmr3"]] <- mu_mr3 + s_mr3 * d_mr3
  rc[["bjb1"]] <- mu_jb1 + s_jb1 * d_jb1
  rc[["bjb3"]] <- mu_jb3 + s_jb3 * d_jb3
  rc[["bc2"]]  <- mu_c2  + s_c2  * d_c2
  rc[["bc3"]]  <- mu_c3  + s_c3  * d_c3
  rc[["bc4"]]  <- mu_c4  + s_c4  * d_c4
  rc[["bcost"]] <- mu_cost + s_cost * tri
  rc[["ec"]]   <- sig_ec * d_ec
  return(rc)
}

apollo_inputs <- apollo_validateInputs()

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P <- list(); V <- list()
  V[["alt1"]] <- b0 + bmr1 * mr1_1 + bmr3 * mr3_1 + bjb1 * jb1_1 + bjb3 * jb3_1 +
                 bc2 * c2_1 + bc3 * c3_1 + bc4 * c4_1 + bcost * cost_1 + ec
  V[["alt2"]] <- b0 + bmr1 * mr1_2 + bmr3 * mr3_2 + bjb1 * jb1_2 + bjb3 * jb3_2 +
                 bc2 * c2_2 + bc3 * c3_2 + bc4 * c4_2 + bcost * cost_2 + ec
  V[["alt3"]] <- 0
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

# ---------------------------------------------------------------------
# 3. WTP for each level (compare with Table 6) and Wald tests
#    WTP = -coefficient / cost coefficient (US$ per household per year)
# ---------------------------------------------------------------------
cat("\n===== WTP by attribute level (compare with Table 6, RPL-EC column) =====\n")
apollo_deltaMethod(model, list(expression = c(
  wtp_b0  = "-mu_b0/mu_cost",
  wtp_mr1 = "-mu_mr1/mu_cost", wtp_mr3 = "-mu_mr3/mu_cost",
  wtp_jb1 = "-mu_jb1/mu_cost", wtp_jb3 = "-mu_jb3/mu_cost",
  wtp_c2  = "-mu_c2/mu_cost",  wtp_c3  = "-mu_c3/mu_cost",  wtp_c4 = "-mu_c4/mu_cost")))

# Gain-loss asymmetry (R2 #4): H0: |WTP(loss)| = WTP(gain)  <=>  WTP(loss) + WTP(gain) = 0.
# A significantly NEGATIVE value means losses weigh more than equal-sized gains.
cat("\n===== Gain-loss asymmetry tests (H0: value = 0) =====\n")
apollo_deltaMethod(model, list(expression = c(
  asym_mr   = "-(mu_mr1 + mu_mr3)/mu_cost",     # WTP(MR -50%) + WTP(MR +50%)
  asym_jobs = "-(mu_jb1 + mu_jb3)/mu_cost")))   # WTP(Jobs -200) + WTP(Jobs +200)

# Differences between carbon levels (R1 #1): H0: no difference in WTP
cat("\n===== Differences between carbon levels (H0: value = 0) =====\n")
apollo_deltaMethod(model, list(expression = c(
  c50_minus_c25  = "-(mu_c3 - mu_c2)/mu_cost",
  c100_minus_c50 = "-(mu_c4 - mu_c3)/mu_cost",
  c100_minus_c25 = "-(mu_c4 - mu_c2)/mu_cost")))
