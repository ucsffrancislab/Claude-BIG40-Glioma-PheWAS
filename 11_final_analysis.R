#!/usr/bin/env Rscript
# =============================================================================
# 11_final_analysis.R
#
# Comprehensive PheWAS analysis, multiple testing correction, sensitivity
# analyses, and publication-ready visualisations.
#
# Runs on EITHER C+T or PRS-CS association results.
#
# Usage:
#   Rscript 11_final_analysis.R <score_dir> <input_dir> <output_dir> [n_cores]
#
# Example:
#   Rscript 11_final_analysis.R ${PWD}/prscs_scores ${PWD}/input ${PWD}/final_results 32
#   Rscript 11_final_analysis.R ${PWD}/ct_output/scores ${PWD}/input ${PWD}/final_results_ct 32
#
# Outputs:
#   <output_dir>/
#     full_results.csv              — All per-cohort results (beta, SE, p, OR, CI)
#     meta_results.csv              — IVW meta-analysis with I², Q, tau²
#     meta_results_corrected.csv    — With M_eff and FDR correction
#     meff_diagnostics.csv          — Eigenvalue summary
#     loocv_results.csv             — Leave-one-cohort-out sensitivity
#     manhattan_<subtype>.png       — PheWAS Manhattan plots
#     forest_top_hits.png           — Forest plots for significant hits
#     heatmap_top_idps.png          — Heatmap of top IDPs × subtypes
#     eigenvalue_scree.png          — Scree plot + M_eff cutoff
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(parallel)
    library(ggplot2)
    library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript 11_final_analysis.R <score_dir> <input_dir> <output_dir> [n_cores]")

SCORE_DIR  <- args[1]
INPUT_DIR  <- args[2]
OUT_DIR    <- args[3]
N_CORES    <- if (length(args) >= 4) as.integer(args[4]) else 32L

COHORTS <- c("cidr", "i370", "onco", "tcga")
PC_COLS <- paste0("PC", 1:8)

OUTCOMES <- list(
    list(col = "case",  label = "all_glioma",  desc = "All glioma vs controls"),
    list(col = "idh",   label = "idh_mutant",  desc = "IDH-mutant vs controls",   val = 1),
    list(col = "idh",   label = "idh_wildtype",desc = "IDH-wildtype vs controls",  val = 0),
    list(col = "pq",    label = "codel",       desc = "1p/19q codel vs controls",  val = 1),
    list(col = "pq",    label = "non_codel",   desc = "1p/19q intact vs controls", val = 0),
    list(col = "tert",  label = "tert_mutant", desc = "TERT-mutant vs controls",   val = 1),
    list(col = "tert",  label = "tert_wildtype",desc= "TERT-wildtype vs controls", val = 0)
)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Load covariates ──────────────────────────────────────────────────────────
cat(sprintf("[%s] Loading covariates ...\n", Sys.time()))
cov_list <- list()
for (cohort in COHORTS) {
    f <- file.path(INPUT_DIR, paste0(cohort, "-covariates.csv"))
    if (!file.exists(f)) { cat(sprintf("  WARNING: %s not found\n", f)); next }
    dt <- fread(f)
    dt[, cohort := cohort]
    if ("exclude" %in% names(dt)) dt <- dt[is.na(exclude) | exclude == 0]
    cov_list[[cohort]] <- dt
    cat(sprintf("  %s: %d samples\n", cohort, nrow(dt)))
}

# ── Discover IDPs ────────────────────────────────────────────────────────────
score_patterns <- c("\\.sscore$", "\\.profile$")
idp_ids <- character(0)
for (pat in score_patterns) {
    files <- list.files(file.path(SCORE_DIR, COHORTS[1]), pattern = pat)
    if (length(files) > 0) {
        idp_ids <- sub(pat, "", files)
        score_ext <- sub("^\\\\", "", sub("\\$$", "", pat))
        break
    }
}
cat(sprintf("[%s] Found %d IDPs (extension: %s)\n", Sys.time(), length(idp_ids), score_ext))

# ── Load IDP metadata if available ───────────────────────────────────────────
idp_meta <- NULL
meta_files <- c("glioma_idps_detail.csv")
for (mf in meta_files) {
    if (file.exists(mf)) { idp_meta <- fread(mf); break }
    mf2 <- file.path(dirname(SCORE_DIR), mf)
    if (file.exists(mf2)) { idp_meta <- fread(mf2); break }
}

# ── Association function (returns beta, SE, p, OR, CI) ───────────────────────
run_association <- function(idp, cohort, outcome) {
    sf <- file.path(SCORE_DIR, cohort, paste0(idp, score_ext))
    if (!file.exists(sf)) return(NULL)
    scores <- tryCatch(fread(sf), error = function(e) NULL)
    if (is.null(scores) || nrow(scores) == 0) return(NULL)

    # Standardize column names
    nms <- names(scores)
    if ("#FID" %in% nms) setnames(scores, "#FID", "FID", skip_absent = TRUE)
    if ("#IID" %in% nms) setnames(scores, "#IID", "IID")
    if (!"IID" %in% names(scores) && "FID" %in% names(scores)) {
        setnames(scores, names(scores)[2], "IID")
    }

    score_col <- names(scores)[ncol(scores)]
    scores[, PGS := as.numeric(get(score_col))]
    scores <- scores[, .(IID, PGS)]

    cov <- cov_list[[cohort]]
    if (is.null(cov)) return(NULL)
    merged <- merge(cov, scores, by = "IID")
    if (nrow(merged) == 0) return(NULL)
    if (sd(merged$PGS, na.rm = TRUE) == 0) return(NULL)
    merged[, PGS := scale(PGS)]

    oc <- outcome
    if (is.null(oc$val)) {
        merged[, Y := as.integer(get(oc$col))]
    } else {
        controls <- merged[case == 0]
        cases_sub <- merged[case == 1 & !is.na(get(oc$col)) & get(oc$col) == oc$val]
        merged <- rbind(controls, cases_sub)
        merged[, Y := as.integer(case)]
    }
    merged <- merged[!is.na(Y)]
    if (length(unique(merged$Y)) < 2) return(NULL)
    n_case <- sum(merged$Y == 1); n_ctrl <- sum(merged$Y == 0)
    if (n_case < 5 || n_ctrl < 5) return(NULL)

    available_pcs <- intersect(PC_COLS, names(merged))
    covars <- "PGS"
    if ("age" %in% names(merged) && any(!is.na(merged$age))) covars <- c(covars, "age")
    if ("sex" %in% names(merged) && any(!is.na(merged$sex))) covars <- c(covars, "sex")
    covars <- c(covars, available_pcs)

    tryCatch({
        fit <- glm(as.formula(paste("Y ~", paste(covars, collapse = "+"))),
                    data = merged, family = binomial())
        coefs <- summary(fit)$coefficients
        if (!"PGS" %in% rownames(coefs)) return(NULL)
        ci <- confint.default(fit)["PGS", ]
        data.table(
            idp = idp, cohort = cohort, outcome = oc$label,
            beta = coefs["PGS", "Estimate"], se = coefs["PGS", "Std. Error"],
            z = coefs["PGS", "z value"], p = coefs["PGS", "Pr(>|z|)"],
            or = exp(coefs["PGS", "Estimate"]),
            ci_lo = exp(ci[1]), ci_hi = exp(ci[2]),
            n_case = n_case, n_ctrl = n_ctrl
        )
    }, error = function(e) NULL, warning = function(w) NULL)
}

# ── Run all associations ─────────────────────────────────────────────────────
cat(sprintf("[%s] Running associations ...\n", Sys.time()))
tasks <- CJ(idp = idp_ids, cohort = names(cov_list), outcome_idx = seq_along(OUTCOMES))
results <- rbindlist(mclapply(seq_len(nrow(tasks)), function(i) {
    run_association(tasks$idp[i], tasks$cohort[i], OUTCOMES[[tasks$outcome_idx[i]]])
}, mc.cores = N_CORES), fill = TRUE)
results <- results[!is.na(beta)]
cat(sprintf("[%s] %d per-cohort results\n", Sys.time(), nrow(results)))
fwrite(results, file.path(OUT_DIR, "full_results.csv"))

# ── Meta-analysis with heterogeneity stats ───────────────────────────────────
cat(sprintf("[%s] Meta-analysis with heterogeneity ...\n", Sys.time()))
meta <- results[, {
    k <- .N
    w <- 1 / se^2
    beta_ivw <- sum(beta * w) / sum(w)
    se_ivw <- sqrt(1 / sum(w))
    z_ivw <- beta_ivw / se_ivw
    p_ivw <- 2 * pnorm(abs(z_ivw), lower.tail = FALSE)
    # Cochran's Q and I²
    Q <- sum(w * (beta - beta_ivw)^2)
    Q_p <- if (k > 1) pchisq(Q, df = k - 1, lower.tail = FALSE) else NA_real_
    I2 <- if (k > 1) max(0, (Q - (k - 1)) / Q * 100) else 0
    tau2 <- if (k > 1) max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w))) else 0
    ci_lo <- exp(beta_ivw - 1.96 * se_ivw)
    ci_hi <- exp(beta_ivw + 1.96 * se_ivw)
    list(beta_meta = beta_ivw, se_meta = se_ivw, z_meta = z_ivw, p_meta = p_ivw,
         or_meta = exp(beta_ivw), ci_lo = ci_lo, ci_hi = ci_hi,
         Q = Q, Q_p = Q_p, I2 = I2, tau2 = tau2,
         n_cohorts = k, total_cases = sum(n_case), total_ctrls = sum(n_ctrl))
}, by = .(idp, outcome)]

# ── M_eff eigenvalue correction ──────────────────────────────────────────────
cat(sprintf("[%s] Computing M_eff (eigenvalue-based correction) ...\n", Sys.time()))

# Build PGS matrix: rows = individuals, cols = IDPs
# Use the first outcome (all_glioma) scores for correlation
pgs_matrix <- NULL
for (cohort in names(cov_list)) {
    for (idp in idp_ids) {
        sf <- file.path(SCORE_DIR, cohort, paste0(idp, score_ext))
        if (!file.exists(sf)) next
        sc <- tryCatch(fread(sf, select = c(ncol(fread(sf, nrows = 0)))), error = function(e) NULL)
        if (is.null(sc)) next
        setnames(sc, 1, idp)
        if (is.null(pgs_matrix)) {
            pgs_matrix <- sc
        } else {
            pgs_matrix <- cbind(pgs_matrix, sc)
        }
        break  # one cohort is enough for correlation structure
    }
}

if (!is.null(pgs_matrix) && ncol(pgs_matrix) > 1) {
    # Standardize
    pgs_mat <- scale(as.matrix(pgs_matrix))
    pgs_mat <- pgs_mat[, apply(pgs_mat, 2, function(x) !any(is.na(x)))]
    
    # Correlation matrix and eigendecomposition
    R <- cor(pgs_mat, use = "pairwise.complete.obs")
    evals <- eigen(R, only.values = TRUE)$values
    
    # Li & Ji M_eff
    M_eff <- sum(ifelse(evals >= 1, 1, 0) + (evals - floor(evals)))
    alpha_eff <- 0.05 / M_eff
    
    cat(sprintf("  N_IDPs = %d, M_eff = %.0f, alpha_eff = %.2e\n",
                ncol(pgs_mat), M_eff, alpha_eff))
    
    # Save eigenvalue diagnostics
    eig_dt <- data.table(
        eigenvalue_rank = seq_along(evals),
        eigenvalue = evals,
        cumvar = cumsum(evals) / sum(evals)
    )
    fwrite(eig_dt, file.path(OUT_DIR, "meff_diagnostics.csv"))
    
    # Eigenvalue scree plot
    fig <- ggplot(eig_dt[1:min(200, nrow(eig_dt))], aes(x = eigenvalue_rank, y = eigenvalue)) +
        geom_point(size = 0.8, alpha = 0.6) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
        annotate("text", x = max(eig_dt$eigenvalue_rank) * 0.7, y = max(evals) * 0.9,
                 label = sprintf("M_eff = %.0f\nalpha = %.2e", M_eff, alpha_eff),
                 hjust = 0, size = 4) +
        labs(title = "PGS Correlation Matrix Eigenvalue Scree Plot",
             x = "Eigenvalue rank", y = "Eigenvalue") +
        theme_minimal()
    ggsave(file.path(OUT_DIR, "eigenvalue_scree.png"), plot = fig, width = 8, height = 5, dpi = 300)
    
} else {
    M_eff <- length(idp_ids)
    alpha_eff <- 0.05 / M_eff
    cat(sprintf("  WARNING: Could not compute PGS correlation matrix. Using naive M_eff = %d\n", M_eff))
}

# Apply corrections
meta[, p_meff := ifelse(p_meta * M_eff > 1, 1, p_meta * M_eff)]
meta[, sig_meff := p_meta < alpha_eff]
meta[, p_fdr := p.adjust(p_meta, method = "fdr"), by = outcome]
meta[, sig_fdr := p_fdr < 0.05]

fwrite(meta, file.path(OUT_DIR, "meta_results_corrected.csv"))

# ── Leave-one-cohort-out ─────────────────────────────────────────────────────
cat(sprintf("[%s] Leave-one-cohort-out sensitivity ...\n", Sys.time()))

# Only for hits with p_meta < 0.01 (save compute)
sig_hits <- unique(meta[p_meta < 0.01, .(idp, outcome)])
if (nrow(sig_hits) > 0) {
    loocv <- rbindlist(lapply(seq_len(nrow(sig_hits)), function(i) {
        hit_idp <- sig_hits$idp[i]
        hit_oc <- sig_hits$outcome[i]
        sub <- results[idp == hit_idp & outcome == hit_oc]
        if (nrow(sub) < 2) return(NULL)
        rbindlist(lapply(unique(sub$cohort), function(drop_cohort) {
            kept <- sub[cohort != drop_cohort]
            if (nrow(kept) == 0) return(NULL)
            w <- 1 / kept$se^2
            b <- sum(kept$beta * w) / sum(w)
            s <- sqrt(1 / sum(w))
            data.table(idp = hit_idp, outcome = hit_oc, dropped = drop_cohort,
                       beta_loocv = b, se_loocv = s,
                       p_loocv = 2 * pnorm(abs(b / s), lower.tail = FALSE),
                       or_loocv = exp(b))
        }))
    }))
    fwrite(loocv, file.path(OUT_DIR, "loocv_results.csv"))
    cat(sprintf("  %d leave-one-out results\n", nrow(loocv)))
} else {
    cat("  No hits at p < 0.01; skipping LOOCV\n")
}

# ── Merge IDP metadata ──────────────────────────────────────────────────────
if (!is.null(idp_meta)) {
    idp_meta[, idp := sprintf("%04d", as.integer(pheno))]
    meta <- merge(meta, idp_meta[, .(idp, idp_short_name, category)],
                  by = "idp", all.x = TRUE)
} else {
    meta[, idp_short_name := idp]
    meta[, category := "unknown"]
}

# ── PheWAS Manhattan Plots ───────────────────────────────────────────────────
cat(sprintf("[%s] Generating Manhattan plots ...\n", Sys.time()))

for (oc_label in unique(meta$outcome)) {
    sub <- meta[outcome == oc_label]
    sub[, idp_num := as.integer(idp)]
    sub <- sub[order(idp_num)]
    sub[, neglogp := -log10(p_meta)]
    
    # Color by category
    cats <- unique(sub$category)
    sub[, cat_idx := as.integer(factor(category, levels = cats))]
    sub[, cat_color := ifelse(cat_idx %% 2 == 0, "A", "B")]
    
    p <- ggplot(sub, aes(x = seq_len(.N), y = neglogp, color = cat_color)) +
        geom_point(size = ifelse(sub$sig_meff, 3, 0.8),
                   alpha = ifelse(sub$neglogp > -log10(0.05), 0.8, 0.3)) +
        geom_hline(yintercept = -log10(alpha_eff), linetype = "solid", color = "red", linewidth = 0.5) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.3) +
        scale_color_manual(values = c("A" = "#1f77b4", "B" = "#ff7f0e"), guide = "none") +
        labs(title = sprintf("PheWAS: %s", oc_label),
             subtitle = sprintf("M_eff = %.0f, threshold = %.2e", M_eff, alpha_eff),
             x = "IDP (ordered by number)", y = expression(-log[10](p))) +
        theme_minimal() +
        theme(axis.text.x = element_blank())
    
    # Label significant points
    sig_sub <- sub[sig_meff == TRUE]
    if (nrow(sig_sub) > 0) {
        p <- p + ggrepel::geom_text_repel(
            data = sig_sub,
            aes(label = idp_short_name),
            size = 2.5, max.overlaps = 20, color = "black"
        )
    }
    
    ggsave(file.path(OUT_DIR, sprintf("manhattan_%s.png", oc_label)),
           plot = p, width = 14, height = 6, dpi = 300)
}

# ── Heatmap of top IDPs ─────────────────────────────────────────────────────
cat(sprintf("[%s] Generating heatmap ...\n", Sys.time()))

# Top 50 IDPs by best p across any outcome
best_per_idp <- meta[, .(best_p = min(p_meta)), by = idp]
top_idps <- best_per_idp[order(best_p)][1:min(50, .N)]$idp

heat_data <- meta[idp %in% top_idps]
heat_data[, signed_neglogp := -log10(p_meta) * sign(beta_meta)]
heat_data[, sig_label := fifelse(sig_meff, "**", fifelse(sig_fdr, "*", ""))]

p_heat <- ggplot(heat_data, aes(x = outcome, y = reorder(idp_short_name, -as.integer(idp)),
                                 fill = signed_neglogp)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sig_label), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, name = "Signed\n-log10(p)") +
    labs(title = "Top 50 IDPs × Glioma Subtypes",
         subtitle = "** = M_eff significant, * = FDR significant",
         x = "Outcome", y = "IDP") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 6),
          axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUT_DIR, "heatmap_top_idps.png"),
       plot = p_heat, width = 10, height = 14, dpi = 300)

# ── Forest plots for significant hits ────────────────────────────────────────
cat(sprintf("[%s] Generating forest plots ...\n", Sys.time()))

sig_idps <- unique(meta[sig_meff == TRUE | sig_fdr == TRUE]$idp)
if (length(sig_idps) == 0) sig_idps <- unique(meta[order(p_meta)][1:min(10, .N)]$idp)

forest_data <- results[idp %in% sig_idps]
forest_data[, ci_lo := exp(beta - 1.96 * se)]
forest_data[, ci_hi := exp(beta + 1.96 * se)]

# Add meta-analysis row
meta_forest <- meta[idp %in% sig_idps, .(idp, outcome, or = or_meta,
                                          ci_lo, ci_hi, cohort = "Meta")]
forest_data_combined <- rbind(
    forest_data[, .(idp, outcome, or = exp(beta), ci_lo, ci_hi, cohort)],
    meta_forest, fill = TRUE
)

for (sig_idp in sig_idps[1:min(10, length(sig_idps))]) {
    fd <- forest_data_combined[idp == sig_idp]
    idp_name <- if (!is.null(idp_meta)) {
        idp_meta[idp == sig_idp]$idp_short_name[1]
    } else sig_idp
    
    p_forest <- ggplot(fd, aes(x = or, y = cohort, xmin = ci_lo, xmax = ci_hi)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_errorbarh(height = 0.2) +
        geom_point(aes(shape = ifelse(cohort == "Meta", "diamond", "circle")), size = 3) +
        scale_shape_identity() +
        facet_wrap(~outcome, nrow = 1) +
        labs(title = sprintf("Forest plot: %s (%s)", sig_idp, idp_name),
             x = "OR (95% CI)", y = "") +
        theme_minimal()
    
    ggsave(file.path(OUT_DIR, sprintf("forest_%s.png", sig_idp)),
           plot = p_forest, width = 14, height = 4, dpi = 300)
}

# ── Summary ──────────────────────────────────────────────────────────────────
cat(sprintf("\n[%s] === RESULTS SUMMARY ===\n", Sys.time()))
cat(sprintf("  M_eff = %.0f (alpha = %.2e)\n\n", M_eff, alpha_eff))

for (oc in unique(meta$outcome)) {
    sub <- meta[outcome == oc]
    n_nom <- sum(sub$p_meta < 0.05)
    n_fdr <- sum(sub$sig_fdr)
    n_meff <- sum(sub$sig_meff)
    best <- sub[which.min(p_meta)]
    cat(sprintf("  %-20s  tested=%d | nominal=%d  FDR=%d  Meff=%d | best: %s (p=%.2e, OR=%.3f)\n",
        oc, nrow(sub), n_nom, n_fdr, n_meff, best$idp, best$p_meta, best$or_meta))
}

cat(sprintf("\n[%s] All outputs saved to %s/\n", Sys.time(), OUT_DIR))
