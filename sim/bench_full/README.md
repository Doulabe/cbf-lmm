# bench_full — Full reproducibility pipeline

End-to-end benchmark and real-data analysis pipeline for the
*Conditional Bayesian Stepwise Refinement for Fine-Mapping in Gaussian
Linear Mixed Models* manuscript.

## Pipeline overview

The pipeline is organised in three layers:

1. **Raw runs** (numeric scripts `01_*.R`–`23_*.R`): each script
   generates simulated data with deterministic seeds, runs the
   relevant methods on each `(cell, replicate)` pair, and writes a
   single RDS per pair to `results/bench_full/<axis>/`.  All scripts
   support `--cores N` for `mclapply` parallelism and are resumable
   (RDS checkpoints).

2. **Aggregation** (`36_aggregate.R`): pools all raw RDS files into
   per-axis summary CSVs and one global aggregate.

3. **Table and figure generators** (`make_*.R`): convert the CSV
   summaries into LaTeX tables and publication figures.

## Anchor configuration

| Parameter | Value |
|---|---|
| K_true | 5 |
| beta_true | (0.8, 0.4, 0.4, 0.2, 0.2) — 1 strong + 2 mod + 2 weak |
| n (anchor) | 1000 |
| m (anchor) | 5000 |
| rho (anchor) | 0.95 (block-AR(1), block_size=10) |
| sigma_g^2 | {0, 0.5} swept on every axis |
| K_max | 10 |
| theta (JointPosterior) | 0.99 |
| N_delta (Gauss-Legendre) | 15 |
| tau^2 (g-prior) | 0.04 |
| B (replicates per cell) | 100 (default), 10–50 for focused/diagnostic scripts |

## Raw-run scripts

| Script | Purpose | Cells × B |
|---|---|---|
| `10_amortized_axis_n.R` | sample size sweep | n ∈ {500, 1000, 3000} × σ² × B=100 |
| `11_amortized_axis_rho.R` | LD strength sweep | ρ ∈ {0.80, 0.95, 0.98} × σ² × B=100 |
| `12_amortized_axis_arch.R` | signal architecture sweep | signal ∈ {weak, med, strong} × σ² × B=100 |
| `13_amortized_axis_m.R` | candidate count sweep | m ∈ {5000, 10000} × σ² × B=100 |
| `90_mouse_full.R` | mouse BMI, chromosome-wise | 20 chromosomes |
| `91_mouse_wholegenome.R` | mouse BMI, whole-genome joint scan | 1 run |
| `24_amortized_prior_sensitivity.R` | half-Cauchy vs inverse-gamma on δ | 10 cells × 2 priors |
| `18_amortized_add_lmm_scan.R` | augment existing RDS with FaST-LMM | all axis dirs |
| `19_amortized_fix_lmm_scan.R` | recover silently-failed FaST-LMM cells | failures only |
| `20_amortized_topk.R` | per-SNP scores for the threshold-independent ranking | cells × B |
| `21_amortized_add_scores.R` | iterative q^Bayes for framework methods | all cells from 10 |
| `22_amortized_topk_remlbf.R` | threshold-independent ranking under REML-BF | cells × B |
| `25_amortized_tau2_sensitivity.R` | τ² ∈ {0.01, 0.04, 0.10, 0.25} | anchor × σ² × B=20 |
| `92_mouse_autosomes.R` | mouse BMI without X chromosome | 1 run |
| `28_amortized_mixtau.R` | mixture-slab variant (negative result) | 5 cells × B=50 |
| `29_amortized_pathavg.R` | path-averaging variant (exploratory) | 2 cells × B=50 |
| `93_mouse_bmi_sex.R` | mouse BMI ± X × ± sex covariate | 2 scans |
| `94_mouse_reml_bf.R` | mouse BMI REML-BF sensitivity scan | 1 run |
| `30_amortized_ablation.R` | variance-ratio × scoring × stopping ablation | anchor × σ² × B |
| `31_rho_V_diagnostic.R` | empirical ρ_V / D_jk distribution under strong LD | anchor |
| `32_bootstrap_ci.R` | bootstrap confidence intervals for the aggregate table | aggregate |
| `23_amortized_threshold_grid.R` | Top-K* recall on the full 16-cell grid | 16 cells |
| `26_eb_tau2_sensitivity.R` | closed-form empirical-Bayes slab (EB-τ²) | 3 cells × B=10 |
| `27_amortized_gamma_sensitivity.R` | eBIC penalty γ ∈ {0, 0.5, 1} | anchor × σ² × B=12 |
| `95_mouse_delta_profile.R` | δ-profile identifiability diagnostic (mouse) | 2 loci |
| `33_shared_kernel_validation.R` | exact vs shared-kernel log-BF (S3 validation) | 4 cells × B=5 |
| `01b–04b_*_remlbf.R` | four-axis benchmark under the REML-BF evaluation | 16 cells |

## Aggregation and reporting

```bash
# Aggregate per-axis raw RDS files into summary CSVs (n ≥ 1000 filter)
Rscript sim/bench_full/36_aggregate.R --filter-n

# Build the main aggregate table (8 methods × 2 polygenic regimes)
Rscript sim/bench_full/70_make_amortized_aggregate_table.R

# Build the runtime table
Rscript sim/bench_full/77_make_amortized_runtime_table.R

# Build per-axis figures (F1, Precision, Recall–FDR)
Rscript sim/bench_full/79_make_amortized_figures.R

# Build axis-stratified supplementary tables
Rscript sim/bench_full/71_make_amortized_supp_tables.R

# Build the inverse-gamma sensitivity table
Rscript sim/bench_full/72_make_amortized_sensitivity_table.R

# Build the τ² sensitivity table
Rscript sim/bench_full/73_make_amortized_tau2_table.R

# Build the threshold-independent comparison table
Rscript sim/bench_full/75_make_amortized_topk_table.R

# Build the mouse BMI 4-way table (chromosomes × sex adjustment)
Rscript sim/bench_full/98_make_mouse_bmi_sex_table.R
```

All generators read from `results/bench_full/` and write to `theory/`.

## One-shot driver

```bash
Rscript sim/bench_full/run_all.R --cores 16 --B 100
```

Runs the four axis sweeps + mouse BMI in sequence, with `mclapply`
parallelism across the cell × replicate grid.  Resumable (existing
RDS checkpoints are skipped).

## Critical environment variables

Always set before running any benchmark script:

```bash
export OMP_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

Without these, Accelerate/OpenBLAS multi-threading collides with
`mclapply` workers and produces ~10–100× slowdowns plus silent worker
failures.  See [`../INSTALL.md`](../../INSTALL.md) for details.

## Total compute budget

| Stage | Approximate CPU-hours |
|---|---|
| 01–04 four-axis benchmark (B=100, 16 cells) | ~400 |
| 01b–04b same axes under REML-BF | ~120 |
| 05–06 mouse BMI (whole-genome + chromosome-wise) | ~2 |
| 07 inverse-gamma sensitivity | ~5 |
| 10–10c threshold-independent + 19 full grid | ~25 |
| 11 τ² + 21 EB-τ² + 22 γ sensitivity | ~15 |
| 12, 15, 15b mouse autosomes ± sex + REML-BF scan | ~4 |
| 13 mixture-τ² benchmark | ~50 |
| 14 path-averaging (exploratory, 2 cells only) | ~10 |
| 16 ablation + 17 ρ_V + 18 bootstrap + 23 δ-profile | ~25 |
| **Total** | **~650 CPU-hours** |

On a 16-core machine with the BLAS thread settings above, this
completes in roughly 40 wall-clock hours.  On a 30-core machine,
roughly 22 wall-clock hours.
