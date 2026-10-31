# Whole-blood expression signature of clinically diagnosed Alzheimer's disease

Code, final model and GEO search for the manuscript *"A whole-blood expression signal of clinically diagnosed Alzheimer's disease transfers across cohorts but is not specific to Alzheimer's disease"* (submitted to *Neurological Sciences*).

A 17-gene whole-blood expression model was developed in GSE140829 and applied without modification to GSE63061 and GSE63060, with random-signature benchmarks, sensitivity analyses and specificity testing against FTLD-spectrum syndromes and MCI (GSE140830). **The signature is not a diagnostic test.**

## Contents

| Path | Contents |
|---|---|
| `code/` | All Python and R analysis scripts, run in order by `run_all.py` |
| `code/R/pipeline_functions.R` | Shared R functions |
| `code/results/frozen_model_coefficients.csv` | Final model: intercept and 17 penalised coefficients (per-SD scale) |
| `code/results/frozen_model_settings.csv` | Classification threshold, penalty and gene-selection rule |
| `code/results/fixed_reference_parameters.csv` | Discovery reference distribution, gene means and SDs for label-free fixed-reference scaling of new arrays |
| `code/results/probe_map_by_cohort.csv`, `signature_probe_identity.csv` | Probe-to-gene mapping in each cohort and identity of signature probe sequences |
| `code/results/checksums_sha256.csv`, `raw_download_checksums.csv` | SHA-256 checksums of raw GEO downloads and of the final model |
| `code/results/python_versions.csv`, `sessionInfo_R.txt` | Software versions |
| `code/data/genesets/Li2014_BTM.gmt` | Blood transcription modules (Li et al., 2014) used in the enrichment analysis |
| `geo_search/` | GEO query, search date and the list of 197 series screened |

## Reproducing the analyses

Requirements: Python 3.11 (numpy, pandas, scipy, statsmodels, matplotlib, requests) and R 4.5 (glmnet, rms, pROC and [EPIC](https://github.com/GfellerLab/EPIC)).

```bash
cd code
python run_all.py
```

The pipeline downloads the public GEO series (GSE140829, GSE140830, GSE63061, GSE63060), performs quality control and regenerates every result and figure. A full run takes about 2 hours on 16 cores. The discovery step must reproduce `results/frozen_model_coefficients.csv` with SHA-256 `6c78e23561798064b9425a8f3f721f14e699b53a05c5417dcd2fc29fee0c31ec`; `run_all.py` prints whether the checksum matches.

## Data

All expression data are public in the Gene Expression Omnibus and are not redistributed here. No participant-level clinical data beyond those deposited in GEO were used.
