# =============================================================================
# BIG40 Glioma PheWAS Pipeline — Complete Guide
# =============================================================================
#
# STATUS (as of Apr 30, 2026):
#   ✅ 08b_prep_target.sh     — VCF → plink binary extraction (done)
#   ✅ 08b_run_ct.sh          — C+T PGS for all 3,935 IDPs × 4 cohorts (done)
#   ✅ 09_ct_association.R    — C+T association testing (done)
#   🔄 08_run_prscs.sh       — PRS-CS on 221 IDPs (cidr/i370/onco ~96%, tcga ~29%)
#   ⏳ 09_compute_pgs.sh     — Score genotypes with PRS-CS weights
#   ⏳ 10_prscs_association.R — PRS-CS association testing
#   ⏳ 11_final_analysis.R    — Full analysis with M_eff, plots, sensitivity
#
# =============================================================================
# STEP-BY-STEP COMMANDS
# =============================================================================

# ── Step 0: Mop up OOM failures from PRS-CS ────────────────────────────────
# Resubmit with fewer CPUs. Crash recovery skips completed IDPs.
for cohort in cidr i370 onco tcga; do
    sbatch --output="logs/prscs-mopup-${cohort}-%j.out" \
           --error="logs/prscs-mopup-${cohort}-%j.err" \
           --mail-user=$(tail -1 ~/.forward) --mail-type=FAIL --time=14-0 \
           --job-name=prscs-mopup-${cohort} \
           --ntasks=1 --cpus-per-task=32 --mem=490G \
           ~/github/ucsffrancislab/Claude-BIG40-Glioma-PheWAS/08_run_prscs.sh \
           ${PWD}/prscs_output "$cohort" ${PWD}/prscs_idp_list.txt
done

# ── Step 1: Score genotypes with PRS-CS weights ────────────────────────────
sbatch --output="logs/pgs-score-%j.out" --error="logs/pgs-score-%j.err" \
       --mail-user=$(tail -1 ~/.forward) --mail-type=FAIL --time=14-0 \
       --job-name=pgs-score \
       --ntasks=1 --cpus-per-task=32 --mem=240G \
       ~/github/ucsffrancislab/Claude-BIG40-Glioma-PheWAS/09_compute_pgs.sh \
       ${PWD}/prscs_output ${PWD}/ct_output/target_bed ${PWD}/prscs_scores \
       ${PWD}/prscs_idp_list.txt

# ── Step 2: Quick PRS-CS association testing ───────────────────────────────
# (Same as C+T association but on PRS-CS scores — just for a quick look)
sbatch --output="logs/prscs-assoc-%j.out" --error="logs/prscs-assoc-%j.err" \
       --mail-user=$(tail -1 ~/.forward) --mail-type=FAIL --time=1-0 \
       --job-name=prscs-assoc --ntasks=1 --cpus-per-task=32 --mem=64G \
       --wrap="module load r; Rscript \
           ~/github/ucsffrancislab/Claude-BIG40-Glioma-PheWAS/10_prscs_association.R \
           ${PWD}/prscs_scores ${PWD}/input 32"

# ── Step 3: Full analysis with M_eff, visualizations, sensitivity ──────────
# Run on PRS-CS scores (primary) and optionally on C+T scores (supplementary)
sbatch --output="logs/final-%j.out" --error="logs/final-%j.err" \
       --mail-user=$(tail -1 ~/.forward) --mail-type=FAIL --time=1-0 \
       --job-name=final --ntasks=1 --cpus-per-task=32 --mem=64G \
       --wrap="module load r; Rscript \
           ~/github/ucsffrancislab/Claude-BIG40-Glioma-PheWAS/11_final_analysis.R \
           ${PWD}/prscs_scores ${PWD}/input ${PWD}/final_results 32"

# Optional: same analysis on C+T scores for supplementary comparison
sbatch --output="logs/final-ct-%j.out" --error="logs/final-ct-%j.err" \
       --mail-user=$(tail -1 ~/.forward) --mail-type=FAIL --time=1-0 \
       --job-name=final-ct --ntasks=1 --cpus-per-task=32 --mem=64G \
       --wrap="module load r; Rscript \
           ~/github/ucsffrancislab/Claude-BIG40-Glioma-PheWAS/11_final_analysis.R \
           ${PWD}/ct_output/scores ${PWD}/input ${PWD}/final_results_ct 32"

# =============================================================================
# WHAT EACH SCRIPT PRODUCES
# =============================================================================
#
# 11_final_analysis.R outputs (in <output_dir>/):
#
#   DATA:
#     full_results.csv              Per-cohort: beta, SE, p, OR, 95% CI
#     meta_results_corrected.csv    IVW meta: + I², Q, tau², M_eff, FDR
#     meff_diagnostics.csv          Eigenvalues for M_eff computation
#     loocv_results.csv             Leave-one-cohort-out sensitivity
#
#   PLOTS:
#     manhattan_<subtype>.png       PheWAS Manhattan (one per subtype)
#     heatmap_top_idps.png          Top 50 IDPs × subtypes heatmap
#     forest_<IDP>.png              Forest plots for significant hits
#     eigenvalue_scree.png          Scree plot with M_eff cutoff
#
# =============================================================================
# PROTOCOL COVERAGE
# =============================================================================
#
# ✅ Phase 1: Data preparation (complete)
# ✅ Phase 2: PGS construction — PRS-CS + C+T
# ✅ Phase 3: Association testing — logistic regression + IVW meta
# ✅ Phase 4: Multiple testing — M_eff eigenvalue + FDR per subtype
# ✅ Phase 5 (partial): Leave-one-cohort-out sensitivity
# ❌ Phase 5 (deferred): Risk-SNP exclusion (Approach A), conditional
#    regression (B), mediation (C), LDSC — these require additional
#    data (glioma GWAS summary stats for PRS) and are best implemented
#    after reviewing initial results to confirm which hits warrant
#    the additional compute.
#
# ✅ Visualisation: Manhattan, forest, heatmap, eigenvalue scree
#
# SUBTYPES TESTED:
#   all_glioma, idh_mutant, idh_wildtype, codel, non_codel,
#   tert_mutant, tert_wildtype
#   (Grade removed per PI recommendation due to inconsistent definitions)
#
# =============================================================================
# SCRIPTS INVENTORY
# =============================================================================
#
#   07c_remap_bim_rsids.py     BIM rsID remapping (done, one-time)
#   08b_prep_target.sh         VCF → plink binary extraction (done)
#   08b_run_ct.sh              C+T scoring all 3,935 IDPs (done)
#   08_run_prscs.sh            PRS-CS posterior weights (running)
#   09_ct_association.R        C+T association testing (done)
#   09_compute_pgs.sh          Score with PRS-CS weights (next)
#   10_prscs_association.R     PRS-CS association testing (next)
#   11_final_analysis.R        Full analysis + plots (next)
#   patch_prscs_numpy2.sh      NumPy 2.x compatibility patch
#   prscs_idp_list.txt         221 IDP list for PRS-CS
#   glioma_idps_detail.csv     IDP annotations (name, category, h²)
#
# =============================================================================
