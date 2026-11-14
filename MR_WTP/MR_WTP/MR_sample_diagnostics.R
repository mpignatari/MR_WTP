# =====================================================================
# Marine reserve DCE - sample diagnostics for JEMECP-D-26-00394
#   Table R1: exclusions by reason and recruitment source (441 -> 371)
#   Table R2: excluded vs retained respondents
#   Table R3: online panel vs mail respondents (analysis sample)
# Base R only. Data: MR_dataset.csv (DigitalCommons)
# NOTE: check value labels in the codebook (Gender, Age, Inc, Educ, Reside, Rec_fish)
# =====================================================================
# setwd("C:/Users/mp1149/Documents/12-2020/Marine Project/JEM submission/reviewer's comments/DATA_Caplan_20260116_ALL/DATA_Caplan_20260116_ALL")
DATAFILE <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "MR_dataset.csv"

raw <- read.csv(DATAFILE)
raw <- raw[order(raw$pid, raw$gid, raw$alt), ]

# respondent-level file + share of status-quo choices (alt 3)
r <- raw[!duplicated(raw$pid), ]
sq <- tapply(raw$Y[raw$alt == 3], raw$pid[raw$alt == 3], mean)
r$sq_share <- as.numeric(sq[as.character(r$pid)])
r$mail     <- as.integer(r$Wave == 4)
r$protest_flag <- r$protest == 1
r$scale_flag   <- r$Excl_sca == 1

# ---------------------------------------------------------------------
# Table R1: exclusions by reason and source
# ---------------------------------------------------------------------
reason <- ifelse(r$protest_flag & r$scale_flag, "Both",
          ifelse(r$protest_flag, "Protest only",
          ifelse(r$scale_flag, "Incomplete attitudinal items only", "Retained")))
reason <- factor(reason, levels = c("Protest only", "Incomplete attitudinal items only",
                                    "Both", "Retained"))
tabR1 <- addmargins(table(reason, Source = ifelse(r$mail == 1, "Mail", "Panel")))
cat("\nTable R1. Exclusions by reason and recruitment source\n"); print(tabR1)
cat("Protest item: rated 6-7 =", sum(r$CE_eval5 %in% 6:7),
    "| missing =", sum(is.na(r$CE_eval5)), "\n")
stopifnot(sum(r$exclall3 == 0) == sum(reason == "Retained"))

# ---------------------------------------------------------------------
# helper: compare two groups
# ---------------------------------------------------------------------
compare <- function(g, lab1, lab0) {
  cat_vars <- c(Gender = "Gender", Age = "Age", Income = "Inc", Education = "Educ",
                Residence = "Reside", Latino = "Latino", `Recreational fishing` = "Rec_fish")
  out <- data.frame(Variable = character(), G1 = character(), G0 = character(),
                    Test = character(), p = numeric(), stringsAsFactors = FALSE)
  for (nm in names(cat_vars)) {
    v <- r[[cat_vars[nm]]]
    ok <- !is.na(v) & !is.na(g)
    p <- suppressWarnings(chisq.test(table(v[ok], g[ok]))$p.value)
    out[nrow(out) + 1, ] <- list(nm,
      sprintf("median cat. %s", median(v[ok & g == 1])),
      sprintf("median cat. %s", median(v[ok & g == 0])), "Chi-square", p)
  }
  if (!all(g == r$mail)) {
    p <- suppressWarnings(chisq.test(table(r$mail, g))$p.value)
    out[nrow(out) + 1, ] <- list("Mail sample (%)",
      sprintf("%.1f", 100 * mean(r$mail[g == 1])), sprintf("%.1f", 100 * mean(r$mail[g == 0])),
      "Chi-square", p)
  }
  p <- wilcox.test(r$sq_share[g == 1], r$sq_share[g == 0], exact = FALSE)$p.value
  out[nrow(out) + 1, ] <- list("Status-quo choices (% of tasks)",
    sprintf("%.1f", 100 * mean(r$sq_share[g == 1])), sprintf("%.1f", 100 * mean(r$sq_share[g == 0])),
    "Mann-Whitney", p)
  p <- suppressWarnings(chisq.test(table(r$sq_share == 1, g))$p.value)
  out[nrow(out) + 1, ] <- list("Always chose status quo (%)",
    sprintf("%.1f", 100 * mean(r$sq_share[g == 1] == 1)),
    sprintf("%.1f", 100 * mean(r$sq_share[g == 0] == 1)), "Chi-square", p)
  names(out)[2:3] <- c(lab1, lab0)
  out$p <- sprintf("%.3f", out$p)
  out
}

# ---------------------------------------------------------------------
# Table R2: excluded vs retained (all 441)
# ---------------------------------------------------------------------
tabR2 <- compare(as.integer(r$exclall3 == 1), "Excluded (70)", "Retained (371)")
cat("\nTable R2. Excluded vs retained respondents\n"); print(tabR2, row.names = FALSE)

# ---------------------------------------------------------------------
# Table R3: panel vs mail (analysis sample, 371)
# ---------------------------------------------------------------------
all_r <- r
r <- r[r$exclall3 == 0, ]
tabR3 <- compare(r$mail, "Mail (59)", "Panel (312)")
cat("\nTable R3. Mail vs online panel respondents (analysis sample)\n"); print(tabR3, row.names = FALSE)
r <- all_r

# category shares for the codebook-labelled version of the tables
cat("\nCategory shares (%), analysis sample, panel vs mail:\n")
for (v in c("Gender", "Age", "Inc", "Educ", "Reside", "Rec_fish")) {
  rr <- r[r$exclall3 == 0, ]
  cat("\n", v, "\n"); print(round(100 * prop.table(table(rr[[v]], ifelse(rr$mail == 1, "Mail", "Panel")), 2), 1))
}
