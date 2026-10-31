
"""run_all.py — reproduce every result, table and figure from public GEO files.
Usage (from this folder):  python run_all.py   [requires Python 3.11 + numpy/pandas/scipy/statsmodels/matplotlib/requests, and Rscript with glmnet/rms/pROC and EPIC (github.com/GfellerLab/EPIC)]
Steps and approximate run times on 16 cores:
  0_download_qc.py            GEO download, QC, harmonisation (~5 min incl. download)
  1_discovery.R               nested CV, final model, 500 bootstraps (~10 min)
  2_validation.R              unchanged-model validation, random signatures, agreement (~3 min)
  0c_fixed_probes.py          probe-sequence identity and probe-level transfer inputs (~2 min)
  3_nomogram_tables.R         nomogram and Tables 4-6 (<1 min)
  0b_scaling_extra_cohorts.py alternative scalings, MCI, GSE140830 (~6 min)
  8_duplicates_ethnicity.py   duplicates GSE140830 vs all pre-QC GSE140829 arrays; AddNeuroMed ethnicity
  5_sensitivity_analyses.R       sensitivity analyses A-F (~25 min)
  7_secondary_analyses.R         secondary analyses: CIs, specificity, benchmarks, supporting statistics (~10 min)
  10_additional_sensitivity_analyses.R  batch-residualised development, subsampling, candidate-pool benchmarks, Brier/loess calibration,
                              age, sex-discordant samples, composition model, EPIC, specificity contrasts, APOE, probe-level transfer (~20 min)
  12_primary_external_validation.R  primary external analysis (fixed-reference scaling, sex-discordant samples retained):
                              performance, benchmarks, sites, sex/age, composition, age/sex-adjusted specificity (~15 min)
  11_enrichment.R             camera / rank tests of Li et al. blood transcription modules (data/genesets/Li2014_BTM.gmt) (<5 min)
  6_site_stratified.py        site-stratified validation; sex-discordant samples by site (<1 min)
  9a_descriptive_tables.py    Tables 1, 1b, signature-gene table, agreement null, checksums, package versions
  4_figures.py                Figures 1-5 and S1
Checks: results/checksums_sha256.csv and results/raw_download_checksums.csv list SHA-256 of raw downloads and of the final model; 1_discovery.R must reproduce
results/frozen_model_coefficients.csv with SHA-256 6c78e23561798064b9425a8f3f721f14e699b53a05c5417dcd2fc29fee0c31ec."""
import subprocess, sys, runpy, hashlib, numpy as np, pandas as pd
def R(script): subprocess.run(["Rscript", "-e", f"source('{script}')"], check=True)
def py(script): runpy.run_path(script, run_name="__main__")
if __name__ == "__main__":
    py("0_download_qc.py"); R("1_discovery.R")
    h = hashlib.sha256(open("results/frozen_model_coefficients.csv", "rb").read()).hexdigest()
    print("final-model checksum match:", h == "6c78e23561798064b9425a8f3f721f14e699b53a05c5417dcd2fc29fee0c31ec")
    py("0c_fixed_probes.py"); R("2_validation.R"); R("3_nomogram_tables.R"); py("0b_scaling_extra_cohorts.py"); py("8_duplicates_ethnicity.py")
    subprocess.run(["Rscript", "-e", "source('5_sensitivity_analyses.R'); source('7_secondary_analyses.R'); source('10_additional_sensitivity_analyses.R'); source('12_primary_external_validation.R'); source('11_enrichment.R')"], check=True)
    py("6_site_stratified.py"); py("9a_descriptive_tables.py"); py("4_figures.py")
