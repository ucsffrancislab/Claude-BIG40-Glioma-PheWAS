# BIG40 Brain Imaging PheWAS — Glioma Subtype Risk
# =============================================================================
#
# Phenome-Wide Association Study testing whether polygenic scores (PGS) for
# 3,935 brain imaging-derived phenotypes (IDPs) from the Oxford BIG40 resource
# predict glioma risk in a subtype-specific manner.
#
# Cohorts: cidr, i370, onco, tcga
# Subtypes: all_glioma, idh_mutant, idh_wildtype, codel, non_codel,
#           tert_mutant, tert_wildtype
#
# Two PGS methods: C+T (fast screen, all 3,935 IDPs) and PRS-CS (Bayesian,
# 221 selected IDPs).
#
# =============================================================================
# PREREQUISITES
# =============================================================================
#
# Cluster modules needed:
#   module load plink    (plink 1.9 — for LD clumping and binary merging)
#   module load plink2   (plink 2.x — for VCF extraction and scoring)
#   module load r        (R with data.table, ggplot2, cowplot, ggrepel)
#
# Input data (expected in ${PWD}/input/):
#   imputed-umich-{cidr,i370,onco,tcga}/chr{1-22}.dose.vcf.gz  (imputed VCFs)
#   {cidr,i370,onco,tcga}-covariates.csv                       (phenotype data)
#
# Covariate CSV columns:
#   IID, dataset, source, age, sex, case, grade, idh, pq, tert, rad, chemo,
#   treated, PC1-PC8, survdays, vstatus, exclude
#
# Reference data (on cluster):
#   BIG40 summary stats:  /francislab/data1/refs/BIG40/prscs_input/
#   Target BIMs (rsID):   /francislab/data1/refs/BIG40/target_bim/
#   PRS-CS LD reference:  ~/.local/ld_ref/ldblk_1kg_eur/
#   PRS-CS software:      ~/.local/PRScs/PRScs.py
#   1000G EUR plink:      /francislab/data1/refs/sources/fileserve.mrcieu.ac.uk/ld/EUR
#
# =============================================================================


# =============================================================================
# STEP 0: Patch PRS-CS for NumPy 2.x compatibility (one-time)
# =============================================================================
# PRS-CS (2021, unmaintained) crashes on NumPy 2.x due to array-to-scalar
# conversion changes. This patch fixes 5 leak points in mcmc_gtb.py.
# Idempotent — safe to rerun. Backs up original as .bak.

bash patch_prscs_numpy2.sh


# =============================================================================
# STEP 1: Prepare target genotypes (one-time, ~1.5 hours)
# =============================================================================
# Extracts EUR-panel SNPs (~8M) from per-chromosome imputed VCFs into plink
# binary filesets. One merged .bed/.bim/.fam per cohort. These are ~100x smaller
# than the raw VCFs and enable fast random-access scoring.
#
# Output: ct_output/target_bed/{cidr,i370,onco,tcga}.{bed,bim,fam}

sbatch --output="logs/ct_prep-%j.out" --error="logs/ct_prep-%j.err" \
       --mail-user=$(tail -1 ~/.forward) --mail-type=FAIL --time=14-0 \
       --job-name=ct_prep \
       --ntasks=1 --cpus-per-task=64 --mem=490G \
       08b_prep_target.sh \
       ${PWD}/input ${PWD}/ct_output


# =============================================================================
# STEP 2: C+T polygenic scoring — all 3,935 IDPs (one-time, ~24 hours)
# =============================================================================
# Clump + Threshold: LD-clump each IDP's summary stats using the 1000G EUR
# panel (plink --clump, r²=0.1, 250kb), then score all 4 cohorts using plink2.
# Clumped SNP lists are cached and write-protected for reuse.
#
# Requires Step 1 (target_bed/) to be complete.
#
# Output: ct_output/scores/{cohort}/{IDP}.profile
#         ct_output/clumped/{IDP}.snps

sbatch --output="logs/ct-%j.out" --error="logs/ct-%j.err" \
       --mail-user=$(tail -1 ~/.forward) --mail-type=FAIL --time=14-0 \
       --job-name=ct \
       --ntasks=1 --cpus-per-task=64 --mem=490G \
       08b_run_ct.sh \
       ${PWD}/ct_output


# =============================================================================
# STEP 3: C+T association testing (~1 minute)
# =============================================================================
# Logistic regression of each C+T PGS against 7 glioma outcomes, per cohort,
# with IVW meta-analysis. Used to select which IDPs go to PRS-CS.
#
# Output: ct_output/association/ct_meta_results.csv
#         ct_output/association/ct_association_summary.csv

sbatch --output="logs/ct-assoc-%j.out" --error="logs/ct-assoc-%j.err" \
       --mail-user=$(tail -1 ~/.forward) --mail-type=FAIL --time=1-0 \
       --job-name=ct-assoc --ntasks=1 --cpus-per-task=64 --mem=64G \
       --wrap="module load r; Rscript \
           09_ct_association.R \
           ${PWD}/ct_output ${PWD}/input 64"


# =============================================================================
# STEP 4: Select IDP subset for PRS-CS
# =============================================================================
# Based on C+T association results: select IDPs with meta-analysis p < 0.01
# in any of the 7 outcomes. This produces the IDP list file.
#
# The current list (prscs_idp_list.txt) contains 221 IDPs.
# A broader list (prscs_idp_list_optionC.txt) with 2,503 IDPs is available
# if compute budget allows.
#
# No command — this was done interactively. To regenerate:
#   python3 -c "
#   import pandas as pd
#   meta = pd.read_csv('ct_output/association/ct_meta_results.csv')
#   hits = sorted(set(str(x).zfill(4) for x in meta[meta['p_meta']<0.01]['idp'].unique()))
#   open('prscs_idp_list.txt','w').write('\n'.join(hits)+'\n')
#   print(f'{len(hits)} IDPs selected')
#   "


# =============================================================================
# STEP 5: PRS-CS posterior weight estimation (~3 days per cohort)
# =============================================================================
# Bayesian MCMC estimation of posterior SNP effect sizes using PRS-CS.
# 500 iterations / 250 burnin for screening (use 1000/500 for final).
# One job per cohort; they queue and run sequentially on the same node.
# Crash recovery: write-protects completed outputs, skips on resubmit.
#
# Output: prscs_output/{cohort}/{IDP}/{IDP}_pst_eff_a1_b0.5_phiauto_chr{1-22}.txt

for cohort in cidr i370 onco tcga; do
    sbatch --output="logs/prscs-${cohort}-%j.out" \
           --error="logs/prscs-${cohort}-%j.err" \
           --mail-user=$(tail -1 ~/.forward) --mail-type=FAIL --time=14-0 \
           --job-name=prscs-${cohort} \
           --ntasks=1 --cpus-per-task=64 --mem=490G \
           08_run_prscs.sh \
           ${PWD}/prscs_output "$cohort" ${PWD}/prscs_idp_list.txt
done

# Mop up OOM failures (resubmit with fewer CPUs after initial run):
for cohort in cidr i370 onco tcga; do
    sbatch --output="logs/prscs-mopup-${cohort}-%j.out" \
           --error="logs/prscs-mopup-${cohort}-%j.err" \
           --mail-user=$(tail -1 ~/.forward) --mail-type=FAIL --time=14-0 \
           --job-name=prscs-mopup-${cohort} \
           --ntasks=1 --cpus-per-task=32 --mem=490G \
           08_run_prscs.sh \
           ${PWD}/prscs_output "$cohort" ${PWD}/prscs_idp_list.txt
done


# =============================================================================
# STEP 6: Score genotypes with PRS-CS weights (~1-2 hours)
# =============================================================================
# Concatenates per-chromosome PRS-CS posterior weights into genome-wide score
# files, translates rsIDs to chr:pos, and scores each cohort using plink2
# against the pre-extracted target binaries.
#
# Output: prscs_scores/{cohort}/{IDP}.sscore

sbatch --output="logs/pgs-score-%j.out" --error="logs/pgs-score-%j.err" \
       --mail-user=$(tail -1 ~/.forward) --mail-type=FAIL --time=14-0 \
       --job-name=pgs-score \
       --ntasks=1 --cpus-per-task=32 --mem=240G \
       09_compute_pgs.sh \
       ${PWD}/prscs_output ${PWD}/ct_output/target_bed ${PWD}/prscs_scores \
       ${PWD}/prscs_idp_list.txt


# =============================================================================
# STEP 7: PRS-CS association testing (~1 minute)
# =============================================================================
# Same logistic regression + IVW meta-analysis as the C+T version (Step 3),
# but using the higher-quality PRS-CS scores.
#
# Output: prscs_scores/association/prscs_meta_results.csv

sbatch --output="logs/prscs-assoc-%j.out" --error="logs/prscs-assoc-%j.err" \
       --mail-user=$(tail -1 ~/.forward) --mail-type=FAIL --time=1-0 \
       --job-name=prscs-assoc --ntasks=1 --cpus-per-task=32 --mem=64G \
       --wrap="module load r; Rscript \
           10_prscs_association.R \
           ${PWD}/prscs_scores ${PWD}/input 32"


# =============================================================================
# STEP 8: Final analysis — M_eff, sensitivity, visualisation (~5-10 minutes)
# =============================================================================
# Comprehensive analysis producing publication-ready outputs:
#   - Eigenvalue-based M_eff multiple testing correction (Li & Ji 2005)
#   - BH-FDR per subtype
#   - Leave-one-cohort-out sensitivity analysis
#   - I², Cochran's Q, tau² heterogeneity statistics
#   - PheWAS Manhattan plots (one per subtype)
#   - Forest plots (for significant hits)
#   - Heatmap (top 50 IDPs × subtypes)
#   - Eigenvalue scree plot
#
# Run on PRS-CS scores (primary analysis):
# Output: final_results/*.csv, final_results/*.png

sbatch --output="logs/final-%j.out" --error="logs/final-%j.err" \
       --mail-user=$(tail -1 ~/.forward) --mail-type=FAIL --time=1-0 \
       --job-name=final --ntasks=1 --cpus-per-task=32 --mem=64G \
       --wrap="module load r; Rscript \
           11_final_analysis.R \
           ${PWD}/prscs_scores ${PWD}/input ${PWD}/final_results 32"

# Optional: same analysis on C+T scores for supplementary comparison:
sbatch --output="logs/final-ct-%j.out" --error="logs/final-ct-%j.err" \
       --mail-user=$(tail -1 ~/.forward) --mail-type=FAIL --time=1-0 \
       --job-name=final-ct --ntasks=1 --cpus-per-task=32 --mem=64G \
       --wrap="module load r; Rscript \
           11_final_analysis.R \
           ${PWD}/ct_output/scores ${PWD}/input ${PWD}/final_results_ct 32"


# =============================================================================
# OPTIONAL: PRS-CS at 1000 iterations for final hits
# =============================================================================
# After reviewing Step 8 results, select the IDPs that survived M_eff or FDR
# correction. Create a refined list and rerun PRS-CS at 1000/500 iterations.
#
# To change iterations, set N_ITER and N_BURNIN before submission:
#   export N_ITER=1000 N_BURNIN=500
# Or edit the defaults in 08_run_prscs.sh (currently 500/250).
#
# Then repeat Steps 6-8 with the refined scores.


# =============================================================================
# DEFERRED ANALYSES (require additional data)
# =============================================================================
#
# These are specified in the protocol but require glioma GWAS summary stats
# or additional computation. Best done after confirming which hits warrant it.
#
# - Glioma risk-SNP exclusion (Approach A): Rebuild PGS excluding ±1Mb around
#   ~25 known glioma risk loci, re-test associations, compute attenuation ratios
#
# - Conditional regression (Approach B): Add glioma PRS as covariate in the
#   association model to test independence from known risk architecture
#
# - Mediation analysis (Approach C): Decompose total/direct/indirect effects
#   of brain PGS on glioma risk through glioma genetic risk
#
# - LDSC genetic correlation: Between top IDP GWAS and glioma subtype GWAS
#
# - Sex-stratified analysis (if sample size permits)


# =============================================================================
# OUTPUT DIRECTORY STRUCTURE
# =============================================================================
#
# ct_output/
#   target_bed/                          Pre-extracted plink binaries (Step 1)
#     {cohort}.{bed,bim,fam}
#   clumped/                             LD-clumped SNP lists (Step 2)
#     {IDP}.snps
#   scores/                              C+T PGS scores (Step 2)
#     {cohort}/{IDP}.profile
#   association/                         C+T association results (Step 3)
#     ct_meta_results.csv
#     ct_association_summary.csv
#     ct_association_results.csv
#   .eur_snps.txt                        EUR panel SNP list (auto-generated)
#   .rsid_to_chrpos.txt                  rsID→chr:pos mapping (auto-generated)
#
# prscs_output/
#   {cohort}/{IDP}/                      PRS-CS posterior weights (Step 5)
#     {IDP}_pst_eff_a1_b0.5_phiauto_chr{1-22}.txt
#   logs/workers/{cohort}/{jobid}/       Per-IDP worker logs
#
# prscs_scores/
#   {cohort}/{IDP}.sscore                PRS-CS PGS scores (Step 6)
#   association/                         PRS-CS association results (Step 7)
#     prscs_meta_results.csv
#
# final_results/                         Publication outputs (Step 8)
#   full_results.csv                     Per-cohort results (beta, SE, p, OR, CI)
#   meta_results_corrected.csv           Meta-analysis + M_eff + FDR
#   meff_diagnostics.csv                 Eigenvalue data
#   loocv_results.csv                    Leave-one-cohort-out
#   manhattan_{subtype}.png              PheWAS Manhattan plots
#   forest_{IDP}.png                     Forest plots
#   heatmap_top_idps.png                 Top 50 IDPs × subtypes
#   eigenvalue_scree.png                 Scree plot


# =============================================================================
# SCRIPTS INVENTORY
# =============================================================================
#
# Script                        Purpose
# ----------------------------  ------------------------------------------------
# patch_prscs_numpy2.sh         NumPy 2.x compatibility patch for PRS-CS
# 07c_remap_bim_rsids.py        BIM SNP ID remapping (chr:pos → rsID)
# 08b_prep_target.sh            VCF → plink binary extraction
# 08b_run_ct.sh                 C+T PGS scoring (all 3,935 IDPs)
# 08_run_prscs.sh               PRS-CS posterior weight estimation
# 09_ct_association.R           C+T association testing
# 09_compute_pgs.sh             Score genotypes with PRS-CS weights
# 10_prscs_association.R        PRS-CS association testing
# 11_final_analysis.R           Full analysis + plots + sensitivity
# prscs_idp_list.txt            221 selected IDPs for PRS-CS
# prscs_idp_list_optionC.txt    2,503 IDPs (broader selection, if needed)
# glioma_idps_detail.csv        IDP metadata (name, category, heritability)


# =============================================================================
# CRASH RECOVERY
# =============================================================================
#
# All SLURM scripts (08b_run_ct.sh, 08_run_prscs.sh, 09_compute_pgs.sh) use
# the same crash-recovery pattern:
#   - Output files are chmod a-w (write-protected) on successful completion
#   - On resubmit, workers check for write-protected outputs and skip
#   - Stale partial outputs (writable) are deleted and recomputed
#   - flock-based locking prevents concurrent workers from double-running
#
# To resume after a crash, timeout, or OOM: resubmit the exact same command.
# No cleanup needed. Completed work is preserved automatically.
