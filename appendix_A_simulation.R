## ============================================================================
## Appendix A - Fully reproducible synthetic demonstration of
## re-randomization-probability weighting (RPW) for population-adjusted
## indirect comparison (PAIC) of binary maintenance outcomes.
##
## Companion code to:
##   Rand K. "Re-randomization-probability weighting for population-adjusted
##   comparisons of maintenance outcomes across trial designs" (manuscript v6;
##   earlier working title: "Harmonizing randomized-withdrawal and
##   treat-through trial designs for population-adjusted indirect comparison
##   of maintenance outcomes using re-randomization-probability weighting").
##
## Copyright (c) 2026 Kim Rand
## SPDX-License-Identifier: MIT
##
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to
## deal in the Software without restriction, including without limitation the
## rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
## sell copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in
## all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
## FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
## IN THE SOFTWARE.
##
## Release note: this file is appendix_A_simulation_v3_kimi.R with only this
## header changed (title, copyright and licence). Code and RNG calls are
## untouched, so every result is identical to the v3 run.
##
## Version 3 (external-review revision): hardened numerics, isolated seeds,
## exact integer comparator counts, general per-patient assignment
## probabilities, and an operating-characteristics section (stratified
## assignment, interval coverage, omitted effect modifier, correlated
## covariates). Terminology: withdrawal-referenced (previously "anchored").
##
## Version 3_kimi (editorial addendum): printouts added so that every number
## quoted in the manuscript is literally produced by this script -- the
## anchoring-decomposition quantities in the demonstration's discussion of the
## withdrawal anchor (new section 6b) and the replicate standard deviation of
## the adjusted withdrawal-referenced estimator in the 2,000-replication Monte
## Carlo. No existing computation or RNG call is touched, so all results are
## identical to the v3 run.
##
## Requirements: base R (>= 4.0); ggplot2 for the figures only.
## The script runs top to bottom from this file's directory and writes
## figures to ./figures/ and a summary table to ./simulation_summary.csv.
##
## Three synthetic trials, loosely inspired by published trials in the atopic
## dermatitis area, but entirely invented:
##   Trial A (IPD available)  : randomized withdrawal; week-16 responders
##                              re-randomized 2:2:1 to high dose : low dose :
##                              withdrawal. Non-responders leave the trial.
##   Trial B (aggregate only) : randomized withdrawal; responders re-randomized
##                              1:2 to withdrawal : high dose. Reports induction
##                              baselines, induction response rate, and
##                              maintenance response among re-randomized
##                              responders -- but NO baselines for responders.
##   Trial C (aggregate only) : treat-through; everyone stays on high dose.
##                              Reports induction baselines, week-16 response,
##                              week-52 response among week-16 responders, and
##                              overall week-52 response (incl. late responders).
## ============================================================================

library(ggplot2)
set.seed(20260807)
dir.create("figures", showWarnings = FALSE)

expit <- function(x) 1 / (1 + exp(-x))
logit <- function(p) log(p / (1 - p))

## ----------------------------------------------------------------------------
## 1. Data-generating model (the simulated "world")
## ----------------------------------------------------------------------------
## Baseline covariates: age (years), male (0/1), sev0 (baseline severity,
## EASI-like scale), prior (prior systemic therapy, 0/1). The trials sample
## from different covariate distributions, so a naive (unadjusted) indirect
## comparison is confounded by population differences.

POPS <- list(
  A = list(age_mean = 38, age_sd = 13, sev_mean = 28, sev_sd = 7, p_male = 0.55, p_prior = 0.40),
  B = list(age_mean = 45, age_sd = 12, sev_mean = 32, sev_sd = 7, p_male = 0.60, p_prior = 0.55),
  C = list(age_mean = 42, age_sd = 13, sev_mean = 30, sev_sd = 7, p_male = 0.50, p_prior = 0.45)
)

draw_covariates <- function(n, pop) {
  data.frame(
    age   = rnorm(n, pop$age_mean, pop$age_sd),
    male  = rbinom(n, 1, pop$p_male),
    sev0  = pmin(72, pmax(16, rnorm(n, pop$sev_mean, pop$sev_sd))),
    prior = rbinom(n, 1, pop$p_prior)
  )
}

## Centred covariates for the outcome models.
zfun <- function(d) cbind(za = (d$age - 40) / 12, zs = (d$sev0 - 30) / 7)

## Induction (week-16) response model. Severity and age are prognostic; the
## severity slope has a drug-specific component (em_sev), i.e. severity is an
## effect modifier -- this is what makes population adjustment matter.
IND <- list(
  A = list(int = -0.45, em_sev = -0.25),
  B = list(int = -0.70, em_sev = -0.10),
  C = list(int = -0.60, em_sev = -0.15)
)
B_AGE <- -0.20; B_SEV <- -0.45; B_PRIOR <- -0.30

p_response <- function(d, drug) {
  z <- zfun(d)
  expit(IND[[drug]]$int + B_AGE * z[, "za"] +
        (B_SEV + IND[[drug]]$em_sev) * z[, "zs"] + B_PRIOR * d$prior)
}

## Maintenance (week-52) model among week-16 responders, by maintenance
## regimen k. Withdrawal is given identical parameters for all drugs, i.e. the
## "equivalent withdrawal properties" transitivity condition for anchored
## comparison holds BY CONSTRUCTION in this simulated world.
MNT_INT <- list(
  A = list(cont_hi = 1.30, cont_lo = 0.75, wd = -0.85),
  B = list(cont_hi = 0.95, wd = -0.85),
  C = list(cont_hi = 1.10)
)
M_AGE <- -0.10; M_SEV <- -0.35

p_maintain <- function(d, drug, k) {
  z <- zfun(d)
  expit(MNT_INT[[drug]][[k]] + M_AGE * z[, "za"] + M_SEV * z[, "zs"])
}

## Late response by week 52 among week-16 non-responders who continue
## treatment (relevant only for the treat-through Trial C).
p_late <- function(d, drug) {
  z <- zfun(d)
  expit(-1.8 + M_SEV * z[, "zs"])
}

## ----------------------------------------------------------------------------
## 2. Generate the three synthetic trials
## ----------------------------------------------------------------------------

## Randomized-withdrawal trial: responders re-randomized to `arms` with
## probabilities `probs`; induction non-responders are not followed further.
make_rr_trial <- function(n, pop, drug, arms, probs) {
  d <- draw_covariates(n, pop)
  d$R <- rbinom(n, 1, p_response(d, drug))
  d$Z <- NA_character_
  resp <- which(d$R == 1)
  d$Z[resp] <- sample(arms, length(resp), replace = TRUE, prob = probs)
  d$M <- 0L
  for (k in arms) {
    sel <- which(d$R == 1 & d$Z == k)
    d$M[sel] <- rbinom(length(sel), 1, p_maintain(d[sel, ], drug, k))
  }
  d$D <- d$R * d$M                      # durable response (week 16 AND week 52)
  d
}

## Treat-through trial: everyone continues on high dose to week 52.
make_tt_trial <- function(n, pop, drug) {
  d <- draw_covariates(n, pop)
  d$R <- rbinom(n, 1, p_response(d, drug))
  d$M <- 0L
  resp <- which(d$R == 1)
  d$M[resp] <- rbinom(length(resp), 1, p_maintain(d[resp, ], drug, "cont_hi"))
  d$D    <- d$R * d$M
  d$late <- 0L
  nonresp <- which(d$R == 0)
  d$late[nonresp] <- rbinom(length(nonresp), 1, p_late(d[nonresp, ], drug))
  d$wk52_any <- pmax(d$D, d$late)       # overall week-52 response, incl. late
  d
}

PI_A <- c(cont_hi = 2/5, cont_lo = 2/5, wd = 1/5)   # Trial A: 2:2:1
PI_B <- c(wd = 1/3, cont_hi = 2/3)                  # Trial B: 1:2

trialA <- make_rr_trial(800, POPS$A, "A", names(PI_A), PI_A)
trialB <- make_rr_trial(700, POPS$B, "B", names(PI_B), PI_B)
trialC <- make_tt_trial(650, POPS$C, "C")

## What the AGGREGATE trials publish (and nothing more): induction-population
## baseline summaries, induction response rate, per-arm maintenance response
## among re-randomized responders (B) or among week-16 responders (C).
publish_rr <- function(d, arms) list(
  n    = nrow(d),
  base = c(age_mean = mean(d$age), age_sd = sd(d$age),
           sev_mean = mean(d$sev0), sev_sd = sd(d$sev0),
           p_male = mean(d$male), p_prior = mean(d$prior)),
  r    = mean(d$R),
  m    = sapply(arms, function(k) mean(d$M[d$R == 1 & d$Z == k])),
  n_k  = sapply(arms, function(k) sum(d$R == 1 & d$Z == k, na.rm = TRUE)),
  s_k  = sapply(arms, function(k) sum(d$M[d$R == 1 & d$Z == k]))   # integer event counts
)
aggB <- publish_rr(trialB, names(PI_B))
aggC <- publish_rr(trialC, NULL)
aggC$n_k <- c(cont_hi = sum(trialC$R == 1))   # treat-through: all responders "continue"
aggC$s_k <- c(cont_hi = sum(trialC$M[trialC$R == 1]))
aggC$m   <- c(cont_hi = mean(trialC$M[trialC$R == 1]))
aggC$wk52_any <- mean(trialC$wk52_any)

## Note the estimand trap: Trial C's OVERALL week-52 response includes late
## responders and is NOT comparable to anything a randomized-withdrawal trial
## can report; the harmonized estimand is durable response r * m.
cat(sprintf("Trial C: overall wk52 response %.3f vs durable response %.3f\n",
            aggC$wk52_any, aggC$r * aggC$m["cont_hi"]))

## ----------------------------------------------------------------------------
## 3. Re-randomization-probability weights and baseline reconstruction
## ----------------------------------------------------------------------------
## Design weight for pseudo-arm k in an IPD randomized-withdrawal trial:
## responders re-randomized to k get weight 1, responders in other arms 0, and
## every induction non-responder gets weight pi_k. In expectation over the
## re-randomization mechanism, every induction patient contributes to pseudo-
## arm k with weight pi_k, so the weighted baseline-covariate distribution
## reproduces that of the full induction population. (This fractional-weight
## form is specific to FIXED re-randomization ratios, as used here. Under
## stratified/unequal assignment do NOT substitute patient-specific pi_ik into
## it; use the inverse-probability form
##   W_ik = (1 - R_i) + R_i * 1(Z_i = k) / pi_ik,   E[W|X,R,S] = 1,
## of which omega_ik = pi_k * W_ik is the fixed-ratio special case.)
rpw_weight <- function(d, k, pi_k) {
  ifelse(d$R == 1, as.numeric(d$Z %in% k), pi_k)
}

wmean <- function(x, w) sum(x * w) / sum(w)
ess   <- function(w) sum(w)^2 / sum(w^2)

## Demonstration: responders are a selected subset (baseline means shifted),
## but each RPW pseudo-arm recovers the induction-population means.
covs <- c("age", "sev0", "male", "prior")
recon <- rbind(
  "Induction population" = colMeans(trialA[covs]),
  "Responders only"      = colMeans(trialA[trialA$R == 1, covs]),
  t(sapply(names(PI_A), function(k) {
    w <- rpw_weight(trialA, k, PI_A[k])
    sapply(trialA[covs], wmean, w = w)
  }))
)
rownames(recon)[3:5] <- paste0("Pseudo-arm ", names(PI_A))
print(round(recon, 3))

## ----------------------------------------------------------------------------
## 4. MAIC weights (Signorovitch method of moments, hand-coded)
## ----------------------------------------------------------------------------
## Weights a_i = exp(X_i alpha) chosen so that the weighted means of the
## matching features equal the comparator's published values. We match the
## means of all four covariates plus the second moment of severity (i.e. its
## SD). Solved by minimizing the convex objective sum(exp(Xc alpha)).
maic_weights <- function(X, target) {
  Xc <- sweep(as.matrix(X), 2, target)
  sds <- apply(Xc, 2, sd)
  if (any(sds <= 0)) stop("zero-variance calibration feature: ", paste(colnames(Xc)[sds <= 0], collapse = ", "))
  Xc   <- sweep(Xc, 2, sds, "/")                # scale for numerical stability
  obj  <- function(al) sum(exp(Xc %*% al))
  grad <- function(al) colSums(Xc * as.vector(exp(Xc %*% al)))
  fit  <- optim(rep(0, ncol(Xc)), obj, grad, method = "BFGS",
                control = list(maxit = 1000, reltol = 1e-12))
  if (fit$convergence != 0) stop("MAIC optimisation did not converge (code ", fit$convergence, ")")
  w <- as.vector(exp(Xc %*% fit$par))
  if (!all(is.finite(w))) stop("non-finite MAIC weights")
  ## exact moment matching check (weighted mean of centred features ~ 0)
  imbal <- max(abs(colSums(Xc * w) / sum(w)))
  attr(w, "max_std_imbalance") <- imbal
  w
}

maic_features <- function(d) cbind(age = d$age, sev0 = d$sev0, sev0sq = d$sev0^2,
                                   male = d$male, prior = d$prior)
## Second-moment target: published sd uses denominator n-1, the weighted second
## moment uses n, so the target is mean^2 + (n-1) sd^2 / n (exact for the
## comparator's empirical second moment).
maic_target <- function(agg) c(
  age    = unname(agg$base["age_mean"]),
  sev0   = unname(agg$base["sev_mean"]),
  sev0sq = unname(agg$base["sev_mean"]^2 + (agg$n - 1) * agg$base["sev_sd"]^2 / agg$n),
  male   = unname(agg$base["p_male"]),
  prior  = unname(agg$base["p_prior"])
)

## ----------------------------------------------------------------------------
## 5. Estimators
## ----------------------------------------------------------------------------
## With calibration weights a (all 1 = no population adjustment):
##   r_hat   : induction response rate, standardized to the target population.
##   q_ht    : durable response under maintenance regimen k, standardized;
##             Horvitz-Thompson form with fixed denominator sum(a) --
##             design-unbiased given a, used for the linear (risk-difference)
##             durable-response estimands.
##   q_rpw   : the pseudo-population ("ratio"/Hajek) form, which conditions on
##             the realized re-randomization counts; consistent for the same
##             quantity and less variable, used for nonlinear (log-odds)
##             contrasts.
##   m_hat   : maintenance conditional on induction response, estimated
##             DIRECTLY as the calibration-weighted success rate among
##             responders actually assigned to arm k. Under the fixed ratios
##             used here, randomization makes these a random subset of all
##             responders, so no pi is needed; under stratified/unequal pi_ik,
##             weight numerator and denominator by 1/pi_ik (Hajek form). The ratio
##             q_ht / r_hat is consistent for the same estimand but carries
##             re-randomization count noise in the numerator that does not
##             cancel against the denominator, inflating its variance.
## pi may be a scalar (fixed ratio) or a per-patient vector of assignment
## probabilities for arm k conditional on response (stratified designs).
estimate_ipd <- function(d, k, pi, a = rep(1, nrow(d))) {
  pi     <- rep(unname(pi), length.out = nrow(d))
  r_hat  <- sum(a * d$R) / sum(a)
  inarm  <- d$R == 1 & d$Z %in% k
  q_ht   <- sum(a[inarm] * d$M[inarm] / pi[inarm]) / sum(a)
  W      <- (1 - d$R) + ifelse(inarm, 1 / pi, 0)      # general IPW form
  q_rpw  <- sum(a * W * d$D) / sum(a * W)
  m_dir  <- sum(a[inarm] * d$M[inarm] / pi[inarm]) / sum(a[inarm] / pi[inarm])
  c(r = r_hat, q_ht = q_ht, q_rpw = q_rpw, m = m_dir, m_ratio = q_ht / r_hat)
}

## Aggregate side: the same estimand family reconstructed from published rates.
estimate_agg <- function(agg, k) {
  q <- agg$r * agg$m[k]
  c(r = agg$r, q = unname(q), m = unname(agg$m[k]))
}

## Anchored contrast through the withdrawal arm, on the log-odds scale of the
## durable-response estimand: Delta_j = logit(q_hi) - logit(q_wd). Nonlinear
## scale, so the conditioned (q_rpw) form is used on the IPD side.
wref_delta <- function(q_hi, q_wd) logit(q_hi) - logit(q_wd)

## Outcome-level pseudo-IPD bootstrap for an aggregate trial: the published
## binary data determine, per patient, one of the cells {induction
## non-responder; responder in arm k with maintenance success; responder in
## arm k with maintenance failure}, so the comparator's outcome-level IPD can
## be reconstructed exactly and resampled (equivalently, one multinomial draw
## over the cell counts). This propagates comparator sampling error, including
## the correlation between the induction response rate and the per-arm
## maintenance rates. The published baseline summaries (the MAIC matching
## targets) are treated as fixed.
resample_agg <- function(agg) {
  succ  <- unname(agg$s_k)                       # exact integer event counts
  fail  <- unname(agg$n_k) - succ
  cells <- c(agg$n - sum(agg$n_k), as.vector(rbind(succ, fail)))
  x <- as.vector(rmultinom(1, agg$n, cells / agg$n))
  K <- length(agg$n_k)
  s <- x[2 * (1:K)]; f <- x[2 * (1:K) + 1]
  agg$r   <- (agg$n - x[1]) / agg$n
  agg$m   <- setNames(s / (s + f), names(agg$n_k))
  agg$n_k <- setNames(s + f, names(agg$n_k))
  agg
}

## ----------------------------------------------------------------------------
## 6. True values of every estimand (possible because the world is synthetic)
## ----------------------------------------------------------------------------
## Expected values computed on a very large covariate draw from the TARGET
## trial's population (outcome probabilities entered analytically, so the only
## error is Monte Carlo error on the covariate integral).
truth_in_pop <- function(pop, drug, k, N = 2e6) {
  d <- draw_covariates(N, pop)
  pr <- p_response(d, drug)
  pm <- p_maintain(d, drug, k)
  r  <- mean(pr); q <- mean(pr * pm)
  c(r = r, q = q, m = q / r)
}
set.seed(1)   # truth integrals independent of trial data
tru_A_hi_inB <- truth_in_pop(POPS$B, "A", "cont_hi")
tru_A_wd_inB <- truth_in_pop(POPS$B, "A", "wd")
tru_B_hi_inB <- truth_in_pop(POPS$B, "B", "cont_hi")
tru_B_wd_inB <- truth_in_pop(POPS$B, "B", "wd")
tru_A_hi_inC <- truth_in_pop(POPS$C, "A", "cont_hi")
tru_C_hi_inC <- truth_in_pop(POPS$C, "C", "cont_hi")

TRUTH <- c(
  dq_AB   = unname(tru_A_hi_inB["q"] - tru_B_hi_inB["q"]),
  dm_AB   = unname(tru_A_hi_inB["m"] - tru_B_hi_inB["m"]),
  wref_AB = unname(wref_delta(tru_A_hi_inB["q"], tru_A_wd_inB["q"]) -
                   wref_delta(tru_B_hi_inB["q"], tru_B_wd_inB["q"])),
  dq_AC   = unname(tru_A_hi_inC["q"] - tru_C_hi_inC["q"])
)

## ----------------------------------------------------------------------------
## 6b. Anchoring-decomposition quantities quoted in the demonstration
## ----------------------------------------------------------------------------
## The "synthetic world" paragraph of the Results quotes standardized
## population quantities, not estimator outputs: (i) the two drugs' true
## conditional maintenance probabilities under withdrawal, standardized to
## Trial B's population, which are equal up to the drugs' different
## responder-selection weights (the withdrawal maintenance parameters are
## identical by construction); (ii) the two induction response rates in the
## same population; (iii) the counterfactual withdrawal-referenced contrast
## with drug B's induction response rate set equal to drug A's, holding every
## maintenance probability fixed; and (iv) the log risk-ratio contrast, on
## which the shared induction-response probability cancels exactly. Printing
## them here makes those quoted numbers reproducible from this script.
cat(sprintf("withdrawal maintenance probabilities, standardized to Trial B's population: drug A %.3f, drug B %.3f\n",
            tru_A_wd_inB["m"], tru_B_wd_inB["m"]))
cat(sprintf("induction response rates, standardized to Trial B's population: drug A %.3f, drug B %.3f\n",
            tru_A_hi_inB["r"], tru_B_hi_inB["r"]))
deltaA <- wref_delta(tru_A_hi_inB["q"], tru_A_wd_inB["q"])
deltaB <- wref_delta(tru_B_hi_inB["q"], tru_B_wd_inB["q"])
## Counterfactual: B's induction response rate set to A's, all maintenance
## behaviour unchanged (q = r * m with r replaced by drug A's r throughout).
deltaB_cf <- wref_delta(tru_A_hi_inB["r"] * tru_B_hi_inB["m"],
                        tru_A_hi_inB["r"] * tru_B_wd_inB["m"])
cat(sprintf("withdrawal-referenced log-OR contrast %.3f; with B's induction rate set equal to A's %.3f; induction-rate component %.3f\n",
            deltaA - deltaB, deltaA - deltaB_cf, (deltaA - deltaB) - (deltaA - deltaB_cf)))
## Same contrast on the log risk-ratio scale: q_jk / q_jwd = m_jk / m_jwd,
## so the (shared-within-regimen) induction-response probability cancels.
cat(sprintf("rate-free log risk-ratio contrast %.3f (direct %.3f)\n",
            log(tru_A_hi_inB["m"] / tru_A_wd_inB["m"]) - log(tru_B_hi_inB["m"] / tru_B_wd_inB["m"]),
            log(tru_A_hi_inB["q"] / tru_A_wd_inB["q"]) - log(tru_B_hi_inB["q"] / tru_B_wd_inB["q"])))

## ----------------------------------------------------------------------------
## 7. Primary analysis: Trial A IPD vs published Trials B and C
## ----------------------------------------------------------------------------
analyse <- function(dA, aggB, aggC) {
  ## MAIC calibration of Trial A's induction population to each target trial
  aB <- maic_weights(maic_features(dA), maic_target(aggB))
  aC <- maic_weights(maic_features(dA), maic_target(aggC))
  one <- rep(1, nrow(dA))

  f <- function(a) list(hi = estimate_ipd(dA, "cont_hi", PI_A["cont_hi"], a),
                        wd = estimate_ipd(dA, "wd",      PI_A["wd"],      a))
  eB <- f(aB); eC <- f(aC); e0 <- f(one)
  bB <- estimate_agg(aggB, "cont_hi"); bBwd <- estimate_agg(aggB, "wd")
  bC <- estimate_agg(aggC, "cont_hi")

  c(## A vs B, unanchored: difference in durable response and in conditional
    ## maintenance, both standardized to Trial B's induction population
    dq_AB_naive = unname(e0$hi["q_ht"] - bB["q"]),
    dq_AB_adj   = unname(eB$hi["q_ht"] - bB["q"]),
    dm_AB_naive = unname(e0$hi["m"]    - bB["m"]),
    dm_AB_adj   = unname(eB$hi["m"]    - bB["m"]),
    dm_AB_ratio_adj = unname(eB$hi["m_ratio"] - bB["m"]),   # q_ht/r form, for CI-width comparison
    ## A vs B, referenced to the withdrawal arms (log-odds scale; conditioned
    ## q_rpw form on the IPD side, see estimate_ipd)
    wref_AB_naive = unname(wref_delta(e0$hi["q_rpw"], e0$wd["q_rpw"]) -
                           wref_delta(bB["q"], bBwd["q"])),
    wref_AB_adj   = unname(wref_delta(eB$hi["q_rpw"], eB$wd["q_rpw"]) -
                           wref_delta(bB["q"], bBwd["q"])),
    ## A vs C (treat-through comparator), unanchored durable response
    dq_AC_naive = unname(e0$hi["q_ht"] - bC["q"]),
    dq_AC_adj   = unname(eC$hi["q_ht"] - bC["q"]),
    ## diagnostics
    ess_B  = ess(aB), ess_B_hi = ess(aB * rpw_weight(dA, "cont_hi", PI_A["cont_hi"])),
    ess_B_wd = ess(aB * rpw_weight(dA, "wd", PI_A["wd"]))
  )
}
primary <- analyse(trialA, aggB, aggC)

## Joint two-sided bootstrap. Each replicate resamples (i) the Trial A IPD at
## the individual level, carrying each patient's joint induction and
## maintenance outcome (non-responders coded as failures) with MAIC weights
## re-estimated inside the replicate, and (ii) the comparators' reconstructed
## outcome-level pseudo-IPD (resample_agg). All correlations within each trial
## -- between induction and maintenance rates, across arms, and between
## numerator and denominator of every ratio -- are thereby preserved, and the
## naive estimators carry uncertainty on the same footing as the adjusted ones.
set.seed(4)   # bootstrap seed isolated from the truth-integration RNG stream
boot <- t(replicate(1000, analyse(trialA[sample(nrow(trialA), replace = TRUE), ],
                                  resample_agg(aggB), resample_agg(aggC))))
nonfinite <- colSums(!is.finite(boot))
if (any(nonfinite > 0)) cat("Non-finite bootstrap replicates (excluded from percentiles):\n",
                            paste(names(nonfinite)[nonfinite > 0], nonfinite[nonfinite > 0], collapse = "; "), "\n")
qfin <- function(x, p) quantile(x[is.finite(x)], probs = p)
ci <- apply(boot[, 1:9], 2, qfin, p = c(0.025, 0.975))

adj_cols   <- c("dq_AB_adj", "dm_AB_adj", "wref_AB_adj", "dq_AC_adj")
naive_cols <- c("dq_AB_naive", "dm_AB_naive", "wref_AB_naive", "dq_AC_naive")
results <- data.frame(
  comparison = c("A vs B durable response, unanchored (risk difference)",
                 "A vs B conditional maintenance, unanchored (risk difference)",
                 "A vs B durable response, withdrawal-referenced (log-OR difference)",
                 "A vs C durable response, unanchored (risk difference)"),
  truth       = TRUTH[c("dq_AB", "dm_AB", "wref_AB", "dq_AC")],
  naive       = primary[naive_cols],
  naive_ci_lo = ci[1, naive_cols],
  naive_ci_hi = ci[2, naive_cols],
  adjusted    = primary[adj_cols],
  ci_lo       = ci[1, adj_cols],
  ci_hi       = ci[2, adj_cols]
)
print(results, row.names = FALSE, digits = 3)
cat(sprintf("dm 95%% PI width: direct (Hajek) %.3f vs q_ht/r ratio form %.3f\n",
            ci[2,"dm_AB_adj"] - ci[1,"dm_AB_adj"],
            ci[2,"dm_AB_ratio_adj"] - ci[1,"dm_AB_ratio_adj"]))
cat(sprintf("max abs standardized imbalance after calibration to B: %.2e\n",
            attr(maic_weights(maic_features(trialA), maic_target(aggB)), "max_std_imbalance")))
cat(sprintf("realized counts, Trial A: responders %d (hi %d, lo %d, wd %d); Trial B: responders %d (wd %d/succ %d, hi %d/succ %d); Trial C: responders %d (succ %d)\n",
            sum(trialA$R), sum(trialA$Z=="cont_hi",na.rm=TRUE), sum(trialA$Z=="cont_lo",na.rm=TRUE), sum(trialA$Z=="wd",na.rm=TRUE),
            sum(trialB$R), aggB$n_k["wd"], aggB$s_k["wd"], aggB$n_k["cont_hi"], aggB$s_k["cont_hi"],
            aggC$n_k["cont_hi"], aggC$s_k["cont_hi"]))
## ESS of the combined calibration-times-design weights per pseudo-arm
## (ESS(w) = (sum w)^2 / sum w^2); this describes covariate information in the
## weighted pseudo-arm, not the outcome information of every Table 3 estimator.
cat(sprintf("ESS (combined a x omega weights): induction MAIC %.0f; pseudo-arm hi %.0f; pseudo-arm wd %.0f (n = %d)\n",
            primary["ess_B"], primary["ess_B_hi"], primary["ess_B_wd"], nrow(trialA)))
write.csv(results, "simulation_summary.csv", row.names = FALSE)

## ----------------------------------------------------------------------------
## 8. Monte Carlo: the adjusted estimators are centred on the truth
## ----------------------------------------------------------------------------
## Trial A is regenerated 2000 times; the comparator aggregates are fixed at
## their population values so that the Monte Carlo isolates the behaviour of
## the IPD-side estimators ("unbiased in expectation" claim).
pop_agg_rr <- function(pop, drug, arms, N = 2e6) {
  d <- draw_covariates(N, pop)
  pr <- p_response(d, drug)
  list(n = N,
       base = c(age_mean = mean(d$age), age_sd = sd(d$age),
                sev_mean = mean(d$sev0), sev_sd = sd(d$sev0),
                p_male = mean(d$male), p_prior = mean(d$prior)),
       r = mean(pr),
       m = sapply(arms, function(k) {
         pm <- p_maintain(d, drug, k); mean(pr * pm) / mean(pr)
       }))
}
set.seed(2)
aggB_pop <- pop_agg_rr(POPS$B, "B", names(PI_B))
aggC_pop <- pop_agg_rr(POPS$C, "C", "cont_hi")

set.seed(3)
mc <- t(replicate(2000, {
  dA <- make_rr_trial(800, POPS$A, "A", names(PI_A), PI_A)
  analyse(dA, aggB_pop, aggC_pop)
}))
mc_cols <- c("dq_AB_naive", "dq_AB_adj", "dm_AB_naive", "dm_AB_adj",
             "wref_AB_naive", "wref_AB_adj", "dq_AC_naive", "dq_AC_adj")
mc_bias <- colMeans(mc[, mc_cols]) -
  rep(TRUTH[c("dq_AB", "dm_AB", "wref_AB", "dq_AC")], each = 2)
cat("Monte Carlo mean bias (naive, adjusted pairs):\n")
print(round(mc_bias, 4))
## Replicate standard deviation quoted in the text for the adjusted
## withdrawal-referenced estimator ("roughly 0.33").
cat(sprintf("replicate SD of wref_AB_adj over 2000 replications: %.3f\n",
            sd(mc[, "wref_AB_adj"])))

## ----------------------------------------------------------------------------
## 9. Figures
## ----------------------------------------------------------------------------
theme_set(theme_minimal(base_size = 11))

## Figure 1: the three maintenance designs -----------------------------------
box <- function(x, y, lab, w = 0.42, h = 0.16) {
  data.frame(xmin = x - w / 2, xmax = x + w / 2,
             ymin = y - h / 2, ymax = y + h / 2, lab = lab, x = x, y = y)
}
mkdesign <- function(trial, boxes, arrows) {
  boxes$trial <- trial; arrows$trial <- trial
  list(boxes = boxes, arrows = arrows)
}
dA_boxes <- rbind(box(0.5, 0.9,  "Induction: Drug A\nn = 800"),
                  box(0.25, 0.55, "Responders"),
                  box(0.78, 0.55, "Non-responders\n(exit trial)"),
                  box(0.14, 0.15, "High dose\n(2)", w = 0.26),
                  box(0.42, 0.15, "Low dose\n(2)", w = 0.26),
                  box(0.70, 0.15, "Withdrawal\n(1)", w = 0.26))
dA_arrows <- data.frame(x = c(0.45, 0.58, 0.20, 0.25, 0.32),
                        y = c(0.82, 0.82, 0.47, 0.47, 0.47),
                        xend = c(0.28, 0.75, 0.14, 0.42, 0.66),
                        yend = c(0.63, 0.63, 0.23, 0.23, 0.23))
dB_boxes <- rbind(box(0.5, 0.9,  "Induction: Drug B\nn = 700"),
                  box(0.25, 0.55, "Responders"),
                  box(0.78, 0.55, "Non-responders\n(exit trial)"),
                  box(0.20, 0.15, "Withdrawal\n(1)", w = 0.30),
                  box(0.55, 0.15, "High dose\n(2)", w = 0.30))
dB_arrows <- data.frame(x = c(0.45, 0.58, 0.20, 0.30),
                        y = c(0.82, 0.82, 0.47, 0.47),
                        xend = c(0.28, 0.75, 0.20, 0.52),
                        yend = c(0.63, 0.63, 0.23, 0.23))
dC_boxes <- rbind(box(0.5, 0.9,  "Induction: Drug C\nn = 650"),
                  box(0.5, 0.15, "All patients continue high dose\n(responders and non-responders)", w = 0.8))
dC_arrows <- data.frame(x = 0.5, y = 0.82, xend = 0.5, yend = 0.24)
des <- list(mkdesign("Trial A (IPD): re-randomization 2:2:1", dA_boxes, dA_arrows),
            mkdesign("Trial B (aggregate): re-randomization 1:2", dB_boxes, dB_arrows),
            mkdesign("Trial C (aggregate): treat-through", dC_boxes, dC_arrows))
des_boxes  <- do.call(rbind, lapply(des, `[[`, "boxes"))
des_arrows <- do.call(rbind, lapply(des, `[[`, "arrows"))
des_boxes$trial  <- factor(des_boxes$trial,  levels = unique(des_boxes$trial))
des_arrows$trial <- factor(des_arrows$trial, levels = levels(des_boxes$trial))
fig1 <- ggplot(des_boxes) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey95", colour = "grey30", linewidth = 0.3) +
  geom_text(aes(x = x, y = y, label = lab), size = 2.7, lineheight = 0.95) +
  geom_segment(data = des_arrows, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = grid::unit(4, "pt")), colour = "grey30",
               linewidth = 0.3) +
  facet_wrap(~trial, nrow = 1) +
  coord_cartesian(xlim = c(0, 1.05), ylim = c(0, 1.02)) +
  theme_void(base_size = 11) + theme(strip.text = element_text(size = 9))
ggsave("figures/fig1_designs.png", fig1, width = 9.5, height = 3, dpi = 600, bg = "white")
ggsave("figures/fig1_designs.pdf", fig1, width = 9.5, height = 3, bg = "white")

## Figure 2: baseline reconstruction by RPW ----------------------------------
sdP <- sapply(trialA[covs], sd)
smd_of <- function(m) (m - colMeans(trialA[covs])) / sdP
smd_df <- do.call(rbind, lapply(rownames(recon)[-1], function(g) {
  data.frame(group = g, covariate = covs, smd = smd_of(recon[g, ]))
}))
smd_df$group <- factor(smd_df$group, levels = rev(rownames(recon)[-1]),
                       labels = rev(c("Responders only", "Pseudo-arm: high dose",
                                      "Pseudo-arm: low dose", "Pseudo-arm: withdrawal")))
smd_df$covariate <- factor(smd_df$covariate, levels = covs,
                           labels = c("Age", "Baseline severity", "Male",
                                      "Prior systemic therapy"))
fig2 <- ggplot(smd_df, aes(x = smd, y = group)) +
  geom_vline(aes(xintercept = 0, linetype = "SMD = 0 (exact agreement)"),
             colour = "grey50") +
  geom_point(aes(shape = "Weighted group mean (SMD vs induction population)"),
             size = 2, colour = "steelblue4") +
  scale_linetype_manual(NULL, values = c("SMD = 0 (exact agreement)" = "solid")) +
  scale_shape_manual(NULL, values = c("Weighted group mean (SMD vs induction population)" = 16)) +
  facet_wrap(~covariate, nrow = 1) +
  labs(x = "Standardized mean difference (SMD) vs full induction population", y = NULL,
       title = "Responders are a selected subset; pseudo-arms approximate the induction population in expectation") +
  theme(legend.position = "bottom") +
  guides(shape = guide_legend(order = 1), linetype = guide_legend(order = 2))
ggsave("figures/fig2_reconstruction.png", fig2, width = 10.5, height = 3, dpi = 600, bg = "white")
ggsave("figures/fig2_reconstruction.pdf", fig2, width = 10.5, height = 3, bg = "white")

## Figure 3: joint bootstrap distributions of all primary estimates ----------
pal <- c("Naive (unadjusted)" = "firebrick", "RPW + MAIC" = "steelblue4")
comp_keys <- c(dq_AB = "A vs B durable response,\nunanchored (risk difference)",
               dm_AB = "A vs B conditional maintenance,\nunanchored (risk difference)",
               wref_AB = "A vs B durable response, withdrawal-\nreferenced (log-OR difference)",
               dq_AC = "A vs C durable response,\nunanchored (risk difference)")
boot_long <- do.call(rbind, lapply(names(comp_keys), function(k) rbind(
  data.frame(comparison = unname(comp_keys[k]), method = "Naive (unadjusted)",
             est = unname(boot[, paste0(k, "_naive")])),
  data.frame(comparison = unname(comp_keys[k]), method = "RPW + MAIC",
             est = unname(boot[, paste0(k, "_adj")])))))
boot_long$comparison <- factor(boot_long$comparison, levels = comp_keys)
## vertical lines: per method, median (solid) and 2.5%/97.5% bounds (dotted)
vlines_of <- function(df) do.call(rbind, lapply(split(df, list(df$comparison, df$method)),
  function(g) data.frame(comparison = g$comparison[1], method = g$method[1],
                         x = unname(quantile(g$est, c(0.5, 0.025, 0.975))),
                         stat = c("Median", "95% percentile interval", "95% percentile interval"))))
vl3 <- vlines_of(boot_long)
tru3 <- data.frame(comparison = factor(comp_keys, levels = comp_keys),
                   truth = TRUTH[names(comp_keys)])
fig3 <- ggplot(boot_long, aes(x = est)) +
  geom_density(aes(fill = method), alpha = 0.45, colour = NA) +
  geom_vline(data = vl3, aes(xintercept = x, colour = method, linetype = stat),
             linewidth = 0.4, show.legend = c(colour = FALSE, linetype = TRUE)) +
  geom_vline(data = tru3, aes(xintercept = truth, linetype = "True value"),
             colour = "black", linewidth = 0.5) +
  facet_wrap(~comparison, scales = "free", nrow = 2) +
  scale_fill_manual(NULL, values = pal) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_linetype_manual(NULL, breaks = c("Median", "95% percentile interval", "True value"),
                        values = c("Median" = "solid", "95% percentile interval" = "dotted",
                                   "True value" = "dashed")) +
  labs(x = "Joint two-sided bootstrap distribution of the estimate (1000 replicates)",
       y = "Density") +
  theme(legend.position = "bottom")
ggsave("figures/fig3_results.png", fig3, width = 9, height = 5.6, dpi = 600, bg = "white")
ggsave("figures/fig3_results.pdf", fig3, width = 9, height = 5.6, bg = "white")

## Figure 4: Monte Carlo sampling distributions ------------------------------
mc_long <- rbind(
  data.frame(est = mc[, "dq_AB_naive"], method = "Naive (unadjusted)",
             comparison = "A vs B durable response (unanchored)"),
  data.frame(est = mc[, "dq_AB_adj"], method = "RPW + MAIC",
             comparison = "A vs B durable response (unanchored)"),
  data.frame(est = mc[, "wref_AB_naive"], method = "Naive (unadjusted)",
             comparison = "A vs B withdrawal-referenced log-OR difference"),
  data.frame(est = mc[, "wref_AB_adj"], method = "RPW + MAIC",
             comparison = "A vs B withdrawal-referenced log-OR difference")
)
tru_mc <- data.frame(comparison = c("A vs B durable response (unanchored)",
                                    "A vs B withdrawal-referenced log-OR difference"),
                     truth = TRUTH[c("dq_AB", "wref_AB")])
vl4 <- vlines_of(mc_long)
fig4 <- ggplot(mc_long, aes(x = est)) +
  geom_density(aes(fill = method), alpha = 0.55, colour = NA) +
  geom_vline(data = vl4, aes(xintercept = x, colour = method, linetype = stat),
             linewidth = 0.4, show.legend = c(colour = FALSE, linetype = TRUE)) +
  geom_vline(data = tru_mc, aes(xintercept = truth, linetype = "True value"),
             colour = "black", linewidth = 0.5) +
  facet_wrap(~comparison, scales = "free") +
  scale_fill_manual(NULL, values = pal) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_linetype_manual(NULL, breaks = c("Median", "95% percentile interval", "True value"),
                        values = c("Median" = "solid", "95% percentile interval" = "dotted",
                                   "True value" = "dashed")) +
  labs(x = "Estimate over 2000 simulated replications",
       y = "Density") +
  theme(legend.position = "bottom")
ggsave("figures/fig4_montecarlo.png", fig4, width = 9, height = 3.6, dpi = 600, bg = "white")
ggsave("figures/fig4_montecarlo.pdf", fig4, width = 9, height = 3.6, bg = "white")

## ----------------------------------------------------------------------------
## 10. Unit checks for the design-weight identities
## ----------------------------------------------------------------------------
## (a) E[W | stratum, R] = 1 under stratified assignment; (b) omega = pi_k * W
## under a fixed ratio.
set.seed(5)
chk <- draw_covariates(2e5, POPS$A); chk$R <- rbinom(nrow(chk), 1, 0.5)
chk$str <- as.integer(chk$sev0 >= 30)
pi_chk <- ifelse(chk$str == 1, 1/3, 2/5)          # arm "hi" probability by stratum
chk$Z <- NA_character_
resp <- which(chk$R == 1)
chk$Z[resp] <- ifelse(runif(length(resp)) < pi_chk[resp], "cont_hi", "other")
W <- (1 - chk$R) + ifelse(chk$R == 1 & chk$Z %in% "cont_hi", 1/pi_chk, 0)
stopifnot(abs(mean(W[chk$str == 0]) - 1) < 0.01, abs(mean(W[chk$str == 1]) - 1) < 0.01)
w_fix <- rpw_weight(trialA, "cont_hi", PI_A["cont_hi"])
W_fix <- (1 - trialA$R) + ifelse(trialA$R == 1 & trialA$Z %in% "cont_hi", 1/unname(PI_A["cont_hi"]), 0)
stopifnot(max(abs(w_fix - unname(PI_A["cont_hi"]) * W_fix)) < 1e-12)
cat("Unit checks passed: E[W|stratum]=1 and omega = pi_k * W (fixed ratio)\n")

## ----------------------------------------------------------------------------
## 11. Operating characteristics and stress scenarios
## ----------------------------------------------------------------------------
## Scenario OC-A: stratified re-randomization. Assignment probabilities vary by
## a baseline-severity stratum (2:2:1 below severity 30; 1:1:1 at or above).
## Correct analysis uses patient-specific pi_ik; the incorrect analysis applies
## the overall marginal share as a single fixed fraction. Outcome models are
## unchanged, so the estimand truths are those already computed.
PI_LOW  <- c(cont_hi = 2/5, cont_lo = 2/5, wd = 1/5)
PI_HIGH <- c(cont_hi = 1/3, cont_lo = 1/3, wd = 1/3)
pi_strat <- function(d, k) ifelse(d$sev0 >= 30, PI_HIGH[k], PI_LOW[k])
make_rr_trial_strat <- function(n, pop, drug) {
  d <- draw_covariates(n, pop)
  d$R <- rbinom(n, 1, p_response(d, drug))
  d$Z <- NA_character_
  for (str in c(FALSE, TRUE)) {
    idx <- which(d$R == 1 & (d$sev0 >= 30) == str)
    pr  <- if (str) PI_HIGH else PI_LOW
    d$Z[idx] <- sample(names(pr), length(idx), replace = TRUE, prob = pr)
  }
  d$M <- 0L
  for (k in names(PI_LOW)) {
    sel <- which(d$R == 1 & d$Z == k)
    d$M[sel] <- rbinom(length(sel), 1, p_maintain(d[sel, ], drug, k))
  }
  d$D <- d$R * d$M
  d
}
set.seed(6)
NS <- 1000
ocA <- t(replicate(NS, {
  dA <- make_rr_trial_strat(800, POPS$A, "A")
  a  <- maic_weights(maic_features(dA), maic_target(aggB_pop))
  pi_ok  <- pi_strat(dA, "cont_hi")
  pi_bad <- mean(pi_strat(dA, "cont_hi")[dA$R == 1])          # marginal share, wrongly applied to all
  bB <- estimate_agg(aggB_pop, "cont_hi")
  c(dq_correct  = unname(estimate_ipd(dA, "cont_hi", pi_ok,  a)["q_ht"] - bB["q"]),
    dq_marginal = unname(estimate_ipd(dA, "cont_hi", pi_bad, a)["q_ht"] - bB["q"]))
}))

## Scenario OC-B: interval coverage of the paired two-trial bootstrap. Each
## replicate draws fresh finite Trials A, B and C, publishes B and C, computes
## the four adjusted contrasts and their percentile intervals (B = 200), and
## records whether each interval covers the truth.
set.seed(7)
NC <- 400; BC <- 200
truths <- unname(TRUTH[c("dq_AB", "dm_AB", "wref_AB", "dq_AC")])
ocB <- t(replicate(NC, {
  dA <- make_rr_trial(800, POPS$A, "A", names(PI_A), PI_A)
  dB <- make_rr_trial(700, POPS$B, "B", names(PI_B), PI_B)
  dC <- make_tt_trial(650, POPS$C, "C")
  gB <- publish_rr(dB, names(PI_B))
  gC <- publish_rr(dC, NULL)
  gC$n_k <- c(cont_hi = sum(dC$R == 1)); gC$s_k <- c(cont_hi = sum(dC$M[dC$R == 1]))
  gC$m <- c(cont_hi = mean(dC$M[dC$R == 1]))
  est <- analyse(dA, gB, gC)
  bt  <- t(replicate(BC, analyse(dA[sample(nrow(dA), replace = TRUE), ],
                                 resample_agg(gB), resample_agg(gC))))
  cols <- c("dq_AB_adj", "dm_AB_adj", "wref_AB_adj", "dq_AC_adj")
  lo <- sapply(cols, function(cc) qfin(bt[, cc], 0.025))
  hi <- sapply(cols, function(cc) qfin(bt[, cc], 0.975))
  se <- sapply(cols, function(cc) sd(bt[is.finite(bt[, cc]), cc]))
  c(est[cols], cover = as.numeric(lo <= truths & truths <= hi), bse = se)
}))

## Scenario OC-C: omitted effect modifier. Calibration matches age, sex and
## prior systemic therapy but omits severity (prognostic and effect-modifying).
feat_red <- function(d) cbind(age = d$age, male = d$male, prior = d$prior)
targ_red <- function(agg) c(age = unname(agg$base["age_mean"]),
                            male = unname(agg$base["p_male"]),
                            prior = unname(agg$base["p_prior"]))
set.seed(8)
ocC <- replicate(NS, {
  dA <- make_rr_trial(800, POPS$A, "A", names(PI_A), PI_A)
  a  <- maic_weights(feat_red(dA), targ_red(aggB_pop))
  unname(estimate_ipd(dA, "cont_hi", PI_A["cont_hi"], a)["q_ht"] - estimate_agg(aggB_pop, "cont_hi")["q"])
})

## Scenario OC-D: correlated covariates with different correlation in the two
## trials (age-severity correlation 0 in Trial A's population, 0.5 in Trial
## B's). Calibration matches the same marginal moments as the main analysis;
## the joint distributions still differ, which marginal moment matching cannot
## remove. Truths are recomputed for this world.
draw_cov_corr <- function(n, pop, rho) {
  z1 <- rnorm(n); z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
  data.frame(age   = pop$age_mean + pop$age_sd * z1,
             male  = rbinom(n, 1, pop$p_male),
             sev0  = pmin(72, pmax(16, pop$sev_mean + pop$sev_sd * z2)),
             prior = rbinom(n, 1, pop$p_prior))
}
set.seed(9)
big <- draw_cov_corr(2e6, POPS$B, 0.5)
q_true_D <- mean(p_response(big, "A") * p_maintain(big, "A", "cont_hi")) -
            mean(p_response(big, "B") * p_maintain(big, "B", "cont_hi"))
aggB_corr <- list(n = 2e6,
                  base = c(age_mean = mean(big$age), age_sd = sd(big$age),
                           sev_mean = mean(big$sev0), sev_sd = sd(big$sev0),
                           p_male = mean(big$male), p_prior = mean(big$prior)),
                  r = mean(p_response(big, "B")),
                  m = c(cont_hi = mean(p_response(big, "B") * p_maintain(big, "B", "cont_hi")) /
                                  mean(p_response(big, "B"))))
rm(big)
set.seed(10)
ocD <- replicate(NS, {
  dA <- make_rr_trial(800, POPS$A, "A", names(PI_A), PI_A)   # rho = 0 in Trial A
  a  <- maic_weights(maic_features(dA), maic_target(aggB_corr))
  unname(estimate_ipd(dA, "cont_hi", PI_A["cont_hi"], a)["q_ht"] - estimate_agg(aggB_corr, "cont_hi")["q"])
})

oc_row <- function(x, truth) c(bias = mean(x) - truth, empSE = sd(x), RMSE = sqrt(mean((x - truth)^2)))
oc <- rbind(
  "OC-A stratified pi, correct W (dq_AB)"          = oc_row(ocA[, "dq_correct"],  TRUTH["dq_AB"]),
  "OC-A stratified pi, marginal fraction (dq_AB)"  = oc_row(ocA[, "dq_marginal"], TRUTH["dq_AB"]),
  "OC-C omitted effect modifier (dq_AB)"           = oc_row(ocC, TRUTH["dq_AB"]),
  "OC-D correlated covariates, marginal matching (dq_AB)" = oc_row(ocD, q_true_D))
covg <- colMeans(ocB[, grep("^cover", colnames(ocB))])
ocB_est <- ocB[, c("dq_AB_adj", "dm_AB_adj", "wref_AB_adj", "dq_AC_adj")]
ocB_bse <- ocB[, grep("^bse", colnames(ocB))]
ocB_sum <- cbind(bias = colMeans(ocB_est) - truths, empSE = apply(ocB_est, 2, sd),
                 RMSE = sqrt(colMeans(sweep(ocB_est, 2, truths)^2)),
                 meanBootSE = colMeans(ocB_bse), coverage = covg)
rownames(ocB_sum) <- c("dq_AB (RD)", "dm_AB (RD)", "wref_AB (log-OR)", "dq_AC (RD)")
cat("\nOC-B: paired two-trial bootstrap, adjusted estimators (", NC, "replicates, B =", BC, "):\n")
print(round(ocB_sum, 4))
cat("\nStress scenarios (", NS, "replicates each):\n")
print(round(oc, 4))
write.csv(ocB_sum, "operating_characteristics_B.csv")
write.csv(oc, "operating_characteristics_stress.csv")

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
cat("Done. Figures in ./figures, summaries in simulation_summary.csv and operating_characteristics_*.csv\n")
