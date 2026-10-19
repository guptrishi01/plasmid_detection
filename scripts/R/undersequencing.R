# undersequencing.R -- side analysis: how likely is a run to be under-sequenced, given its
# sequencing effort and what kind of sample it is?
#   Rscript undersequencing.R <nonpareil_all.tsv> <samples.tsv> <site_groups> <threshold>
#                             <sweep> <pilot.tsv> <out_dir>
#   site_groups: "name:SITE,SITE name:SITE ..." (params.stats.site_groups)
#   sweep:       "0.4 0.5 0.6 0.7"
#   pilot.tsv:   per pilot run: raw-read and host-depleted Nonpareil coverage, assembled %
#
# Outcome: Nonpareil modelled coverage C (share of the community the reads have sampled).
# Under-sequenced = C < threshold. Two models, same predictors:
#   logistic  P(C < t) ~ log2(effort) + site_group + system        (binary outcome)
#   firth     the same, Firth-penalised (logistf): finite estimates when a site group is
#             entirely on one side of the threshold (separation), which plain logistic can't
#   beta      C        ~ log2(effort) + site_group + system        (C kept continuous; logit
#             link, so exp(coef) is the change in the odds C/(1-C); P(C < t) from the fitted
#             beta distribution)
# Diversity (Nd) is NOT a predictor: it comes from the same curve as C, so it would be
# circular. Effort = bp in the reads Nonpareil saw (raw R1).
# Compared by leave-one-out cross-validation at every sweep threshold: AUC, Brier score,
# log loss of the predicted P(C < t). Outputs in <out_dir>.
suppressMessages({
  library(betareg)
  library(logistf)
  library(ggplot2)
})

a <- commandArgs(trailingOnly = TRUE)
np <- read.delim(a[1])
sheet <- read.delim(a[2])
groups_spec <- strsplit(a[3], " ")[[1]]
threshold <- as.numeric(a[4])
sweep <- as.numeric(strsplit(a[5], " ")[[1]])
pilot <- read.delim(a[6])
out <- a[7]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ---- data --------------------------------------------------------------------------
site_group <- c()
for (g in groups_spec) {
  p <- strsplit(g, ":")[[1]]
  for (s in strsplit(p[2], ",")[[1]]) site_group[s] <- p[1]
}
d <- merge(np, sheet[, c("run", "location", "site", "timepoint", "spots")], by = "run")
stopifnot(nrow(d) == nrow(np))
d$site_group <- unname(site_group[d$site])
if (any(is.na(d$site_group))) {
  stop("sites with no group: ", paste(unique(d$site[is.na(d$site_group)]), collapse = ", "))
}
# Reference level = the largest group, so the odds ratios read against the most data
d$site_group <- relevel(factor(d$site_group), names(which.max(table(d$site_group))))
d$location <- factor(d$location)
d$effort_gbp <- d$effort_bp / 1e9
d$log2_effort <- log2(d$effort_gbp)
eps <- 1e-4
d$coverage_b <- pmin(pmax(d$coverage, eps), 1 - eps)   # betareg needs 0 < C < 1

# Location has 4 levels: the two creek systems plus the two remote controls, which are
# already their own site_group. "system" = the creek system; controls are put on the
# reference level so the control site_group term carries them (no collinear column).
creeks <- sort(setdiff(unique(as.character(d$location[d$site_group != "control"])), NA))
d$system <- factor(ifelse(d$site_group == "control", creeks[1], as.character(d$location)), levels = creeks)
form_rhs <- if (length(creeks) > 1) "log2_effort + site_group + system" else "log2_effort + site_group"
write.table(d, file.path(out, "undersequencing_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

fit_logit <- function(df, t) {
  df$under <- as.integer(df$coverage < t)
  suppressWarnings(glm(as.formula(paste("under ~", form_rhs)), data = df, family = binomial))
}
fit_firth <- function(df, t) {
  df$under <- as.integer(df$coverage < t)
  logistf(as.formula(paste("under ~", form_rhs)), data = df)
}
predict_firth <- function(m, newdata) {
  X <- model.matrix(delete.response(terms(m$formula)), newdata,
                    xlev = lapply(m$data[, sapply(m$data, is.factor), drop = FALSE], levels))
  plogis(drop(X %*% coef(m)))
}
fit_beta <- function(df) betareg(as.formula(paste("coverage_b ~", form_rhs)), data = df, link = "logit")
p_under_beta <- function(m, newdata, t) {
  mu <- predict(m, newdata, type = "response")
  phi <- predict(m, newdata, type = "precision")
  pbeta(t, mu * phi, (1 - mu) * phi)
}
auc <- function(y, p) {
  if (length(unique(y)) < 2) return(NA)
  r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
logloss <- function(y, p) {
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

# ---- leave-one-out comparison at every threshold -----------------------------------
cv <- do.call(rbind, lapply(sweep, function(t) {
  y <- as.integer(d$coverage < t)
  if (length(unique(y)) < 2) {
    return(data.frame(threshold = t, n_under = sum(y), n = length(y), model = c("logistic", "firth", "beta"),
                      auc = NA, brier = NA, log_loss = NA, note = "one class only: not fittable"))
  }
  pl <- pf <- pb <- numeric(nrow(d))
  for (i in seq_len(nrow(d))) {
    tr <- d[-i, ]; te <- d[i, , drop = FALSE]
    # a factor level seen only in the held-out run can't be predicted: fall back to its
    # training base rate
    ok <- all(sapply(c("site_group", "system"), function(v) te[[v]] %in% tr[[v]]))
    if (!ok) { pl[i] <- pf[i] <- pb[i] <- mean(tr$coverage < t); next }
    pl[i] <- predict(fit_logit(tr, t), te, type = "response")
    pf[i] <- predict_firth(fit_firth(tr, t), te)
    pb[i] <- p_under_beta(fit_beta(tr), te, t)
  }
  sep <- tryCatch(any(fitted(fit_logit(d, t)) < 1e-6 | fitted(fit_logit(d, t)) > 1 - 1e-6),
                  error = function(e) TRUE)
  data.frame(threshold = t, n_under = sum(y), n = length(y), model = c("logistic", "firth", "beta"),
             auc = c(auc(y, pl), auc(y, pf), auc(y, pb)),
             brier = c(mean((pl - y)^2), mean((pf - y)^2), mean((pb - y)^2)),
             log_loss = c(logloss(y, pl), logloss(y, pf), logloss(y, pb)),
             note = c(if (sep) "separation: estimates unreliable" else "", "", ""))
}))
write.table(cv, file.path(out, "model_comparison_loocv.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- full fits at the main threshold ------------------------------------------------
y <- as.integer(d$coverage < threshold)
coefs <- list()
if (length(unique(y)) == 2) {
  ml <- fit_logit(d, threshold)
  ci <- suppressMessages(tryCatch(confint.default(ml), error = function(e) NULL))
  s <- summary(ml)$coefficients
  coefs$logistic <- data.frame(model = "logistic", term = rownames(s),
    odds_ratio = exp(s[, 1]), ci_low = exp(ci[, 1]), ci_high = exp(ci[, 2]), p = s[, 4],
    meaning = "odds of being under-sequenced (C < threshold)")
  mf <- fit_firth(d, threshold)
  coefs$firth <- data.frame(model = "firth", term = names(coef(mf)),
    odds_ratio = exp(coef(mf)), ci_low = exp(mf$ci.lower), ci_high = exp(mf$ci.upper), p = mf$prob,
    meaning = "odds of being under-sequenced (C < threshold), Firth-penalised")
  sink(file.path(out, "logistic_summary.txt")); print(summary(ml))
  cat("\nconverged:", ml$converged, " fitted probabilities at 0/1:",
      sum(fitted(ml) < 1e-6 | fitted(ml) > 1 - 1e-6), "\n")
  cat("\n==== Firth-penalised logistic (logistf) ====\n"); print(summary(mf)); sink()
}
mb <- fit_beta(d)
sb <- summary(mb)$coefficients$mean
cib <- confint(mb)[rownames(sb), , drop = FALSE]
coefs$beta <- data.frame(model = "beta", term = rownames(sb),
  odds_ratio = exp(sb[, 1]), ci_low = exp(cib[, 1]), ci_high = exp(cib[, 2]), p = sb[, 4],
  meaning = "odds of coverage C/(1-C); < 1 = less covered")
sink(file.path(out, "beta_summary.txt")); print(summary(mb)); sink()
coef_tab <- do.call(rbind, coefs)
write.table(coef_tab, file.path(out, "odds_ratios.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Per-run predicted probability of under-sequencing (beta model, main threshold) and the
# effort each run would need to reach the threshold, from Nonpareil's own curve
d$p_under_beta <- p_under_beta(mb, d, threshold)
d$under <- y
write.table(d[order(-d$p_under_beta), c("run", "location", "site", "site_group", "timepoint",
            "effort_gbp", "coverage", "under", "p_under_beta", "diversity", "effort_bp_for_95")],
            file.path(out, "per_run_predictions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- plots --------------------------------------------------------------------------
g1 <- ggplot(d, aes(effort_gbp, coverage, colour = site_group)) +
  geom_hline(yintercept = threshold, linetype = 2, colour = "grey40") +
  geom_point(size = 2.4) + scale_x_log10() +
  labs(x = "Sequencing effort (Gbp, raw R1, log scale)", y = "Nonpareil coverage (C)",
       colour = "Site group", title = "Coverage vs effort, all 92 runs",
       subtitle = sprintf("Dashed line: under-sequencing threshold C = %.2f", threshold)) +
  theme_bw()
ggsave(file.path(out, "coverage_vs_effort.png"), g1, width = 8, height = 5, dpi = 150)

grid <- expand.grid(effort_gbp = exp(seq(log(min(d$effort_gbp) / 2), log(max(d$effort_gbp) * 4), length.out = 100)),
                    site_group = levels(d$site_group))
grid$system <- factor(creeks[1], levels = creeks)
grid$site_group <- factor(grid$site_group, levels = levels(d$site_group))
grid$log2_effort <- log2(grid$effort_gbp)
grid$p <- p_under_beta(mb, grid, threshold)
g2 <- ggplot(grid, aes(effort_gbp, p, colour = site_group)) + geom_line(linewidth = 1) +
  geom_vline(xintercept = range(d$effort_gbp), linetype = 3, colour = "grey50") +
  scale_x_log10() + ylim(0, 1) +
  labs(x = "Sequencing effort (Gbp, log scale)", y = sprintf("P(coverage < %.2f)", threshold),
       colour = "Site group", title = "Predicted probability of under-sequencing (beta regression)",
       subtitle = sprintf("System: %s. Dotted: observed effort range", creeks[1])) +
  theme_bw()
ggsave(file.path(out, "p_undersequenced_vs_effort.png"), g2, width = 8, height = 5, dpi = 150)

# ---- pilot: raw-read vs host-depleted coverage, and coverage vs assembled fraction ----
write.table(pilot, file.path(out, "pilot_calibration.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("runs:", nrow(d), " under-sequenced at C <", threshold, ":", sum(y), "\n")
print(cv[, c("threshold", "n_under", "model", "auc", "brier", "log_loss", "note")], row.names = FALSE)
print(coef_tab[, c("model", "term", "odds_ratio", "ci_low", "ci_high", "p")], row.names = FALSE)
