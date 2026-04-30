#!/usr/bin/env Rscript
# =============================================================================
# 10_prscs_association.R
#
# Association testing of PRS-CS polygenic scores against glioma outcomes.
# Tests each IDP's PGS against all available phenotype outcomes across
# all 4 cohorts, with inverse-variance weighted meta-analysis.
#
# Usage:
#   Rscript 10_prscs_association.R <ct_output_dir> <input_dir> [n_cores]
#
# Example:
#   Rscript 10_prscs_association.R ${PWD}/ct_output ${PWD}/input 32
#
# Inputs:
#   <ct_output_dir>/scores/<cohort>/<IDP>.profile  (from 09_compute_pgs.sh)
#   <input_dir>/<cohort>-covariates.csv             (IID, case, age, sex, PCs, ...)
#
# Outputs:
#   <ct_output_dir>/association/prscs_association_results.csv   (all results)
#   <ct_output_dir>/association/prscs_association_summary.csv   (top hits)
#   <ct_output_dir>/association/prscs_meta_results.csv          (meta-analysis)
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript 10_prscs_association.R <ct_output_dir> <input_dir> [n_cores]")
}

SCORE_BASE    <- args[1]
INPUT_DIR <- args[2]
N_CORES   <- if (length(args) >= 3) as.integer(args[3]) else 32L

COHORTS  <- c("cidr", "i370", "onco", "tcga")
PC_COLS  <- paste0("PC", 1:8)

# Outcomes: each is a list(col, type, description)
# Binary outcomes: logistic regression
# For subtypes: compare subtype-positive cases vs controls
OUTCOMES <- list(
    list(col = "case",    label = "all_glioma",      desc = "All glioma vs controls"),
    list(col = "idh",     label = "idh_mutant",      desc = "IDH-mutant cases vs controls",    val = 1),
    list(col = "idh",     label = "idh_wildtype",    desc = "IDH-wildtype cases vs controls",  val = 0),
    list(col = "pq",      label = "codel",           desc = "1p/19q codel cases vs controls",  val = 1),
    list(col = "pq",      label = "non_codel",       desc = "1p/19q intact cases vs controls", val = 0),
    list(col = "tert",    label = "tert_mutant",     desc = "TERT-mutant cases vs controls",   val = 1),
    list(col = "tert",    label = "tert_wildtype",   desc = "TERT-wildtype cases vs controls", val = 0)
)

cat(sprintf("[%s] Loading covariate files ...\n", Sys.time()))

# ── Load covariates ──────────────────────────────────────────────────────────
load_covariates <- function(cohort) {
    f <- file.path(INPUT_DIR, paste0(cohort, "-covariates.csv"))
    if (!file.exists(f)) {
        cat(sprintf("  WARNING: %s not found, skipping cohort\n", f))
        return(NULL)
    }
    dt <- fread(f)
    dt[, cohort := cohort]
    # Apply exclusion filter
    if ("exclude" %in% names(dt)) {
        n_before <- nrow(dt)
        dt <- dt[is.na(exclude) | exclude == 0]
        cat(sprintf("  %s: %d samples (%d excluded)\n", cohort, nrow(dt), n_before - nrow(dt)))
    } else {
        cat(sprintf("  %s: %d samples\n", cohort, nrow(dt)))
    }
    dt
}

cov_list <- lapply(COHORTS, load_covariates)
names(cov_list) <- COHORTS
cov_list <- cov_list[!sapply(cov_list, is.null)]

# ── Discover available IDPs ──────────────────────────────────────────────────
score_dir <- SCORE_BASE
first_cohort <- COHORTS[1]
idp_files <- list.files(file.path(score_dir, first_cohort), pattern = "\\.sscore$")
idp_ids <- sub("\\.sscore$", "", idp_files)
cat(sprintf("[%s] Found %d IDPs to test\n", Sys.time(), length(idp_ids)))

# ── Association function ─────────────────────────────────────────────────────
run_association <- function(idp, cohort, outcome) {
    # Load score file
    score_file <- file.path(score_dir, cohort, paste0(idp, ".sscore"))
    if (!file.exists(score_file)) return(NULL)

    scores <- tryCatch(fread(score_file), error = function(e) NULL)
    if (is.null(scores) || nrow(scores) == 0) return(NULL)

    # Standardize column names
    # plink2 .sscore: #FID IID ALLELE_CT NAMED_ALLELE_DOSAGE_SUM SCORE1_SUM
    # or:             #IID ALLELE_CT NAMED_ALLELE_DOSAGE_SUM SCORE1_SUM
    if ("#FID" %in% names(scores)) {
        setnames(scores, c("#FID", "IID"), c("FID", "IID"), skip_absent = TRUE)
    } else if ("#IID" %in% names(scores)) {
        setnames(scores, "#IID", "IID")
    }

    # Get the score column (last column)
    score_col <- names(scores)[ncol(scores)]
    scores[, PGS := as.numeric(get(score_col))]
    scores <- scores[, .(IID, PGS)]

    # Merge with covariates
    cov <- cov_list[[cohort]]
    if (is.null(cov)) return(NULL)
    merged <- merge(cov, scores, by = "IID")
    if (nrow(merged) == 0) return(NULL)

    # Standardize PGS to mean=0, sd=1
    if (sd(merged$PGS, na.rm = TRUE) == 0) return(NULL)
    merged[, PGS := scale(PGS)]

    # Build outcome variable
    oc <- outcome
    if (is.null(oc$val)) {
        # Direct binary (case column)
        merged[, Y := as.integer(get(oc$col))]
    } else {
        # Subtype: keep controls (case==0) and cases with specific value
        controls <- merged[case == 0]
        cases_sub <- merged[case == 1 & !is.na(get(oc$col)) & get(oc$col) == oc$val]
        merged <- rbind(controls, cases_sub)
        merged[, Y := as.integer(case)]
    }

    merged <- merged[!is.na(Y)]
    if (length(unique(merged$Y)) < 2) return(NULL)

    # Check we have enough samples
    n_case <- sum(merged$Y == 1)
    n_ctrl <- sum(merged$Y == 0)
    if (n_case < 5 || n_ctrl < 5) return(NULL)

    # Build formula with available PCs
    available_pcs <- intersect(PC_COLS, names(merged))
    covars <- c("PGS")
    if ("age" %in% names(merged) && any(!is.na(merged$age))) covars <- c(covars, "age")
    if ("sex" %in% names(merged) && any(!is.na(merged$sex))) covars <- c(covars, "sex")
    covars <- c(covars, available_pcs)

    formula_str <- paste("Y ~", paste(covars, collapse = " + "))

    # Fit logistic regression
    result <- tryCatch({
        fit <- glm(as.formula(formula_str), data = merged, family = binomial())
        coefs <- summary(fit)$coefficients
        if (!"PGS" %in% rownames(coefs)) return(NULL)
        data.table(
            idp      = idp,
            cohort   = cohort,
            outcome  = oc$label,
            beta     = coefs["PGS", "Estimate"],
            se       = coefs["PGS", "Std. Error"],
            z        = coefs["PGS", "z value"],
            p        = coefs["PGS", "Pr(>|z|)"],
            n_case   = n_case,
            n_ctrl   = n_ctrl,
            or       = exp(coefs["PGS", "Estimate"])
        )
    }, error = function(e) NULL, warning = function(w) NULL)

    result
}

# ── Run all associations in parallel ─────────────────────────────────────────
cat(sprintf("[%s] Running associations (%d IDPs × %d cohorts × %d outcomes = %d tests, %d cores) ...\n",
    Sys.time(), length(idp_ids), length(cov_list), length(OUTCOMES),
    length(idp_ids) * length(cov_list) * length(OUTCOMES), N_CORES))

tasks <- CJ(idp = idp_ids, cohort = names(cov_list), outcome_idx = seq_along(OUTCOMES))

results <- mclapply(seq_len(nrow(tasks)), function(i) {
    run_association(tasks$idp[i], tasks$cohort[i], OUTCOMES[[tasks$outcome_idx[i]]])
}, mc.cores = N_CORES)

results <- rbindlist(results[!sapply(results, is.null)])
cat(sprintf("[%s] %d association results\n", Sys.time(), nrow(results)))

# ── Meta-analysis (inverse-variance weighted) ────────────────────────────────
cat(sprintf("[%s] Running IVW meta-analysis across cohorts ...\n", Sys.time()))

meta <- results[, {
    w <- 1 / se^2
    beta_ivw <- sum(beta * w) / sum(w)
    se_ivw   <- sqrt(1 / sum(w))
    z_ivw    <- beta_ivw / se_ivw
    p_ivw    <- 2 * pnorm(abs(z_ivw), lower.tail = FALSE)
    n_cohorts <- .N
    total_cases <- sum(n_case)
    total_ctrls <- sum(n_ctrl)
    list(
        beta_meta   = beta_ivw,
        se_meta     = se_ivw,
        z_meta      = z_ivw,
        p_meta      = p_ivw,
        or_meta     = exp(beta_ivw),
        n_cohorts   = n_cohorts,
        total_cases = total_cases,
        total_ctrls = total_ctrls
    )
}, by = .(idp, outcome)]

# Multiple testing correction
meta[, p_bonferroni := p.adjust(p_meta, method = "bonferroni"), by = outcome]
meta[, p_fdr := p.adjust(p_meta, method = "fdr"), by = outcome]

# ── Save results ─────────────────────────────────────────────────────────────
out_dir <- file.path(SCORE_BASE, "association")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fwrite(results, file.path(out_dir, "prscs_association_results.csv"))
fwrite(meta,    file.path(out_dir, "prscs_meta_results.csv"))

# Summary: top hits per outcome
top <- meta[order(p_meta)][, head(.SD, 20), by = outcome]
fwrite(top, file.path(out_dir, "prscs_association_summary.csv"))

cat(sprintf("\n[%s] === RESULTS SUMMARY ===\n", Sys.time()))
for (oc in unique(meta$outcome)) {
    sub <- meta[outcome == oc]
    n_nom  <- sum(sub$p_meta < 0.05)
    n_fdr  <- sum(sub$p_fdr < 0.05)
    n_bonf <- sum(sub$p_bonferroni < 0.05)
    best   <- sub[which.min(p_meta)]
    cat(sprintf("  %-20s  %d tested | nominal=%d  FDR=%d  Bonf=%d | best: IDP %s (p=%.2e, OR=%.3f)\n",
        oc, nrow(sub), n_nom, n_fdr, n_bonf, best$idp, best$p_meta, best$or_meta))
}

cat(sprintf("\n[%s] Results saved to %s/\n", Sys.time(), out_dir))
cat(sprintf("  prscs_association_results.csv  (%d rows — per-cohort results)\n", nrow(results)))
cat(sprintf("  prscs_meta_results.csv         (%d rows — meta-analysis)\n", nrow(meta)))
cat(sprintf("  prscs_association_summary.csv  (top 20 per outcome)\n"))
