# CBF-LMM: Conditional Bayes-Factor Selection with an Adaptive Polygenic Background

**CBF-LMM** is a conditional Bayes-factor stepwise procedure for **refining**
candidate regions under an **adaptive polygenic background**: previously
selected SNPs enter as fixed effects and, for each candidate under
evaluation, the remaining markers *other than that candidate* define its
linear mixed-model polygenic background (model kernel `K_jk`), rebuilt along
the selection path.  The primary implementation (`R/CBF_LMM_exact.R`)
evaluates this **exact candidate-specific model at shared-decomposition
cost**: one pool eigendecomposition per step plus Sherman–Morrison rank-one
downdates deliver every `V_jk^{-1}` and `|V_jk|` in closed form (validated
against a brute-force per-candidate eigendecomposition to 2e-13).
The candidate enters as one additional fixed effect; the residual variance is integrated out analytically, the
variance-component ratio $\delta$ is **numerically marginalized** by
one-dimensional Gauss–Legendre quadrature, and model growth is controlled by
a **profile-ML extended BIC** on the full accepted model (the empty set is a
possible return).  A per-candidate **REML plug-in** for $\delta$, a
**Score-LMM ranking ablation** (`rank_by = "score"`), a posterior-score
stopping rule, and the **amortized pool-kernel variant**
(`MS_L_LMM_stepwise_fast()`) are available as sensitivity options. 

> Reference: K. N. Doulabe and L. Lakhal-Chaieb, *Conditional Bayes-Factor
> Selection with an Adaptive Polygenic Background*, 2026 (submitted).

The method is positioned as a **second-stage refinement tool** downstream of a
genome-wide screen or a marginal Bayesian fine-mapper (e.g. SuSiE), not as a
replacement for them: it returns a compact, parsimony-controlled set of
candidates.

## Paper-to-code name mapping

| Paper | Code | Status |
|---|---|---|
| **CBF-LMM** (primary: exact `K_jk` + profile-ML eBIC) | `CBF_LMM_stepwise_exact()` in `R/CBF_LMM_exact.R` (defaults: `delta_eval = "marginal"`) | primary |
| Per-candidate REML plug-in for $\delta$ | `CBF_LMM_stepwise_exact(delta_eval = "reml")` | sensitivity |
| Conditional Score-LMM ranking ablation | `CBF_LMM_stepwise_exact(rank_by = "score")` | ablation |
| Amortized pool-kernel variant (`G_k` for all candidates) | `CBF_LMM_stepwise()` = `MS_L_LMM_stepwise_fast()` | variant (Appendix B) |
| Posterior-score stopping | `criterion = "JointPosterior"` (amortized variant) | sensitivity |


## Highlights

- **Closed-form conditional Bayes factor** at fixed $\delta$, valid for any
  step $k \ge 1$ under a point-normal $g$-prior on the candidate effect, a
  Jeffreys prior on $\sigma^2$, and weakly-informative priors on the nuisance
  parameters.  Strictly generalizes the Wakefield approximate Bayes factor
  (recovered as $\delta \to 0$ at $k = 1$).
- **Primary operating point (CBF-LMM)**: $\delta$-marginalized evaluation
  (Gauss–Legendre quadrature over a half-Cauchy prior on $\sqrt\delta$) with
  eBIC stopping ($\gamma = 1$, parsimony-oriented).
- **Sensitivity options**: a fast plug-in **REML** evaluation of $\delta$
  (justified by an $O_p(n^{-1/2})$ asymptotic-equivalence result and in close
  empirical agreement at fine-mapping sample sizes) and a posterior-score
  stopping rule (`criterion = "JointPosterior"`, recall-oriented).
- **Covariate extension**: the conditioning block extends to arbitrary fixed
  covariates (sex, batch, family, …) without modification to the Bayes-factor
  algebra. This extension and the exploratory **joint Schur-complement** score
  (more stable under within-block LD saturation) are discussed as future
  directions in the manuscript.
- **Comparators** (same $X$, $y$): BayesB (BGLR; 2,000 iterations / 400 burn-in)
  and BayesR (hibayes; 10,000 iterations / 2,000 burn-in, per the chain-length
  control of Appendix D), sparse Bayesian whole-genome regression mixtures; SuSiE (credible-set); and a
  single-kernel LMM association scan in the spirit of FaST-LMM (rrBLUP).
- **Sensitivity suite**: slab variance $\tau^2$, a closed-form empirical-Bayes
  slab (EB-$\tau^2$), the eBIC penalty $\gamma$, the prior on $\delta$
  (half-Cauchy vs inverse-gamma), and a $\delta$-profile identifiability
  diagnostic on the mouse panel.
- **Reproducibility**: full benchmark pipeline with seed-deterministic
  cell-level checkpointing.

## Quick start

```r
# Load the framework
source("R/LMM_core.R")
source("R/LMM_stepwise_fast.R")
source("R/LMM_reml_bf.R")        # needed for the delta_eval = "reml" sensitivity option
source("R/CBF_LMM.R")            # paper-facing alias

# Simulate a small dataset: n=500, m=300, K_true=3 active SNPs
n <- 500L; m <- 300L; K_true <- 3L
set.seed(42)
X <- matrix(rnorm(n * m), n, m)
X <- scale(X)
truth <- c(50L, 120L, 200L)
beta_true <- c(0.8, 0.4, 0.4)
y <- as.numeric(X[, truth] %*% beta_true) + rnorm(n)

# Primary procedure (CBF-LMM): marginalised delta + eBIC stopping — the defaults.
res <- CBF_LMM_stepwise(y, X, tau2 = 0.04, K_max = 10L, n_nodes = 15L)
cat("Selected indices:", res$indices, "\n")
cat("Truth:           ", truth, "\n")
cat("K_hat:", res$K_hat, " (true K_true =", K_true, ")\n")

# Sensitivity option: plug-in REML evaluation of delta
res_reml <- CBF_LMM_stepwise(y, X, tau2 = 0.04, K_max = 10L,
                             delta_eval = "reml")
```

See [`examples/quickstart.R`](examples/quickstart.R) for a complete runnable
example. 

### Choosing the slab variance $\tau^2$ (empirical-Bayes option)

The default $\tau^2 = 0.04$ corresponds to a prior standard deviation of $0.2$
residual-scale units for a standardized candidate effect; in applications it
can be anchored to the anticipated standardized effect size.  Because the
conditional Bayes factor is available in closed form in $\tau^2$, a per-step
**empirical-Bayes estimate** of the slab along the selection path is provided:

```bash
Rscript sim/bench_full/26_eb_tau2_sensitivity.R
```

The script ranks candidates at the default, maximizes the closed-form Bayes
factor in $\tau^2$ for the leading candidate at each step (`tau2_mode = "eb"`
in the in-script stepwise), reports the distribution of the per-step estimates,
and verifies that the selected sets match the fixed default on the anchor
design (Appendix D of the manuscript).

## Repository layout

```
.
├── README.md                  This file
├── INSTALL.md                 Detailed installation / dependency notes
├── LICENSE                    MIT
├── CITATION.cff               Citation metadata
├── R/                         Core framework implementations
│   ├── LMM_core.R                Closed-form Bayes-factor primitives
│   ├── LMM_stepwise.R            Per-candidate reference implementation
│   ├── LMM_stepwise_fast.R       Batched version (delta_eval = "marginal"/"reml")
│   ├── LMM_reml_bf.R             Plug-in REML-BF evaluation (delta_eval = "reml")
│   ├── LMM_stepwise_mixtau.R     Mixture-slab variant (sensitivity)
│   ├── LMM_stepwise_pathavg.R    Path-averaged variant (exploratory)
│   └── CBF_LMM.R                 Paper-facing alias CBF_LMM_stepwise()
├── sim/                       Simulation and benchmark scripts
│   ├── bench_full/               Full reproducibility pipeline
│   │   ├── 00_config.R              Common configuration + comparator runners
│   │   ├── 01-04_*.R                Four-axis benchmark (n, ρ, signal, m)
│   │   ├── 01b-04b_*_remlbf.R       Same four axes under the REML-BF evaluation
│   │   ├── 05-06_*.R                Mouse BMI (chromosome-wise + whole-genome)
│   │   ├── 24_amortized_prior_sensitivity.R   Half-Cauchy vs inverse-gamma on δ
│   │   ├── 18_amortized_add_lmm_scan.R         Add the single-kernel LMM scan to existing RDS
│   │   ├── 20_amortized_topk.R, 10c_*  Threshold-independent ranking
│   │   ├── 25_amortized_tau2_sensitivity.R    τ² sensitivity
│   │   ├── 12_mouse_autosomes_*.R   Autosomes-only mouse BMI scan
│   │   ├── 28_amortized_mixtau.R    Mixture-slab benchmark
│   │   ├── 29_amortized_pathavg.R   Path-averaging benchmark (exploratory)
│   │   ├── 93_mouse_bmi_sex.R  Mouse BMI with sex covariate
│   │   ├── 94_mouse_reml_bf.R  REML-BF real-data sensitivity scan
│   │   ├── 30_amortized_ablation.R            Variance-ratio × scoring × stopping ablation
│   │   ├── 31_rho_V_diagnostic.R    Empirical ρ_V / D_jk distribution under LD
│   │   ├── 32_bootstrap_ci.R        Bootstrap confidence intervals
│   │   ├── 23_amortized_threshold_grid.R Top-K* recall on the full 16-cell grid
│   │   ├── 26_eb_tau2_sensitivity.R Closed-form empirical-Bayes slab (EB-τ²)
│   │   ├── 27_amortized_gamma_sensitivity.R   eBIC penalty γ sensitivity
│   │   ├── 95_mouse_delta_profile.R δ-profile identifiability diagnostic
│   │   ├── 33_shared_kernel_validation.R  Shared pool G_k vs model kernel K_jk (per step)
│   │   ├── 37_semisynth_1000g.R    Semi-synthetic benchmark on real 1000G panels
│   │   ├── 36_aggregate.R           Pool & summarise raw RDS files
│   │   ├── make_*.R                 LaTeX table / figure generators
│   │   ├── 66_make_cbf_figures.R       CBF-LMM manuscript figures (5-method roster)
│   │   ├── 67_make_cbf_rho_table.R     CBF-LMM LD-axis supplement table
│   │   ├── 38_geuvadis_loci_scan.R       GEUVADIS human cis-eQTL illustration
│   │   ├── 65_make_geuvadis_table.R    GEUVADIS results table
│   │   ├── 34_path_kernel_validation.R  Full-path validation of the shared-pool amortisation
│   │   ├── 40_exact_axes.R         Exact-K_jk rerun of the controlled benchmark (same seeds)
│   │   ├── 42_exact_reml_reference.R   REML-plug-in arm on the anchor (exact engine)
│   │   ├── 43_exact_topk.R          Top-K* ranking with the exact engine (forced K*)
│   │   ├── 44_exact_semisynth.R     Semi-synthetic 1000G rerun (exact engine)
│   │   ├── 45_exact_geuvadis.R      GEUVADIS rerun (exact engine)
│   │   ├── 96_mouse_exact.R         Mouse autosomes rerun (exact engine)
│   │   ├── 60_make_exact_tables.R   Manuscript tables from the exact campaign
│   │   ├── 61_make_summaries_figures.R  Summary CSVs + manuscript figures
│   │   ├── 62_make_pooled_tables.R       THE pooling rule (16 unique cells, anchor once)
│   │   ├── 46_exact_tau2_gamma.R  Exact slab/eBIC-penalty sweeps on the anchor
│   │   ├── 47_exact_score_ablation.R     Bayes-factor vs Q^2 ranking ablation
│   │   └── run_all.R                Single-command driver (axes 01–06)
│   └── validate_lmm_bayes.R      Six-check internal validation (V1–V6)
├── examples/                  Minimal demos
│   ├── quickstart.R              20-line single-dataset demo
│   ├── simulation_example.R      One anchor cell with full pipeline
│   ├── real_data_example.R       Minimal mouse-BMI analysis
│   └── delta_profile_diagnostic.R  δ-profile QC plot (REML-BF vs MBF rule)
└── tests/                     Sanity checks (R)
```

## Reproducing the benchmark

The simulation grid has **16 unique cells** spanning four axes — sample size
$n \in \{500, 1000, 3000\}$, block-AR(1) correlation $\rho \in \{0.80, 0.95,
0.98\}$, signal architecture $\in \{$weak, medium, strong$\}$, and number of
SNPs $m \in \{5000, 10000\}$ — each crossed with the homogeneous-noise and
polygenic regimes $\sigma_g^2 \in \{0, 0.5\}$.  Everything reproduces from
([`sim/bench_full/`](sim/bench_full/)):

```bash
# Stage 0: minimal reproducible example (no external data; a few minutes on a laptop)
Rscript examples/minimal_example.R

# Stage 1x: amortized pool-kernel campaign (variant used by the exploratory arms of Appendix D)
Rscript sim/bench_full/run_all.R                                   # 10-13: four design axes, in parallel
Rscript sim/bench_full/14_amortized_axis_n_remlbf.R                # 14-17: REML-BF arms on the same axes
Rscript sim/bench_full/15_amortized_axis_rho_remlbf.R
Rscript sim/bench_full/16_amortized_axis_arch_remlbf.R
Rscript sim/bench_full/17_amortized_axis_m_remlbf.R
Rscript sim/bench_full/18_amortized_add_lmm_scan.R                 # single-kernel LMM scan comparator

# Stage 2x: ranking and hyperparameter sensitivities (amortized variant)
Rscript sim/bench_full/20_amortized_topk.R                         # Top-K* recall on the 16-cell grid
Rscript sim/bench_full/23_amortized_threshold_grid.R
Rscript sim/bench_full/24_amortized_prior_sensitivity.R            # half-Cauchy vs inverse-gamma on delta
Rscript sim/bench_full/25_amortized_tau2_sensitivity.R             # slab variance tau^2
Rscript sim/bench_full/26_eb_tau2_sensitivity.R                    # empirical-Bayes slab update (Appendix E)
Rscript sim/bench_full/27_amortized_gamma_sensitivity.R            # eBIC penalty gamma

# Stage 3x: ablations, validations, real-data preparation
Rscript sim/bench_full/30_amortized_ablation.R
Rscript sim/bench_full/32_bootstrap_ci.R
Rscript sim/bench_full/33_shared_kernel_validation.R               # per-step ranking agreement (Appendix B)
Rscript sim/bench_full/34_path_kernel_validation.R                 # full-path prefix/stopping agreement (Appendix B)
Rscript sim/bench_full/36_aggregate.R --filter-n
Rscript sim/bench_full/37_semisynth_1000g.R --B 200                # requires the two locus panels, see below
Rscript sim/bench_full/38_geuvadis_loci_scan.R                     # GEUVADIS loci, see the data section below
Rscript sim/bench_full/39_geuvadis_susie_cs.R                      # SuSiE credible sets at the GEUVADIS loci

# Stage 4x: exact K_jk campaign (primary results of the paper)
Rscript sim/bench_full/40_exact_axes.R --cores 5                   # controlled evaluation, same seeds as stage 1x
Rscript sim/bench_full/41_exact_aggregate.R
Rscript sim/bench_full/42_exact_reml_reference.R                   # REML plug-in arm (reference configuration)
Rscript sim/bench_full/43_exact_topk.R --cores 5                   # Top-K* ranking
Rscript sim/bench_full/44_exact_semisynth.R --cores 5              # semi-synthetic 1000G
Rscript sim/bench_full/45_exact_geuvadis.R                         # GEUVADIS
Rscript sim/bench_full/46_exact_tau2_gamma.R                       # slab and eBIC-penalty sweeps (Appendix D)
Rscript sim/bench_full/47_exact_score_ablation.R                   # Bayes factor vs Q^2 ranking ablation (Appendix D)
Rscript sim/bench_full/48_exact_reml_rho098.R                      # same-node ablation arm at rho = 0.98

# Stage 5x: robustness and sensitivity analyses (Appendices C-D)
Rscript sim/bench_full/50_semisynth_matching.R --cores 5 --B 200   # one-to-one matching metrics (Table 4)
Rscript sim/bench_full/51_pip_thresholds.R     --cores 5 --B 100   # SuSiE/BayesR PIP-threshold sensitivity (Appendix D)
Rscript sim/bench_full/52_mcmc_stability.R     --cores 5 --B 10    # MCMC chain-length control (Appendix D)
Rscript sim/bench_full/53_correlated_causals.R --cores 5 --B 100   # same-block correlated-causal stress test (Appendix C)
Rscript sim/bench_full/54_delta_bounds.R       --cores 5 --B 100   # delta-grid truncation check (Appendix B)
bash sim/bench_full/58_driver_bayesr_long.sh                       # 55-57: BayesR at 10,000/2,000 across every analysis
                                                                   #        (in place; short-chain rows kept as BayesR_2k)

# Stage 6x: manuscript tables and figures
Rscript sim/bench_full/62_make_pooled_tables.R                     # pooled tables (paper pooling rule)
Rscript sim/bench_full/60_make_exact_tables.R                      # remaining manuscript tables
Rscript sim/bench_full/61_make_summaries_figures.R                 # summary CSVs and figures
Rscript sim/bench_full/63_make_matching_tables.R                   # Table 4 and by-locus table
Rscript sim/bench_full/64_make_appendix_tables.R                   # Appendix C/D tables
Rscript sim/bench_full/65_make_geuvadis_table.R
Rscript sim/bench_full/66_make_cbf_figures.R
Rscript sim/bench_full/67_make_cbf_rho_table.R

# Stage 7x: tables and figures of the amortized variant (not in the paper)
# Stage 9x: mouse BMI analyses (not reported in the paper)
```

### Semi-synthetic 1000G benchmark

Script 37 (and its exact counterpart, 44) evaluates the same methods on **real** genotypes (real LD) rather
than block-AR(1) simulations.  The two locus panels are **not redistributed
here**; they are built from the public 1000 Genomes phase-3 release.  Each panel
is an `.rds` holding a list with

| field | content |
|---|---|
| `G` | `n × m` genotype matrix, integer dosages `0/1/2` |
| `R` | `m × m` LD matrix (`cor(G)`) |
| `pos` | physical positions, length `m` |
| `snpvar` | per-SNP prior variance (unused by this script) |

The published run uses the EUR subset (`n = 503`, MAF ≥ 0.05) at a chromosome-1
locus (`m = 1,493`) and a chromosome-6 locus (`m = 2,143`).  Place them at
`data/1000g/locus_chr1_1000g.rds` and `data/1000g/locus_chr6_1000g.rds`, or point
the script elsewhere:

```bash
Rscript sim/bench_full/37_semisynth_1000g.R \
    --locus_chr1 /path/to/chr1.rds --locus_chr6 /path/to/chr6.rds --B 200
```

Recovery is scored in the LD-aware sense standard for real-LD fine-mapping (a
causal counts as recovered, and a selection as a true positive, at
`r² ≥ 0.25`), because at `n = 503` a causal variant routinely has a perfect-LD
twin and exact-index recovery is not identifiable.

### GEUVADIS cis-eQTL illustration

Step 9 fine-maps eight strong cis-eQTL genes in the GEUVADIS LCL panel
(n = 358 individuals shared with the 1000 Genomes phase-3 EUR reference),
using the published EUR373 best-association list as a concordance benchmark.
All inputs are public; genotypes are streamed remotely per locus by
`bcftools` (required on PATH):

```bash
mkdir -p data/geuvadis
B=http://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEUV/E-GEUV-1/analysis_results
curl -o data/geuvadis/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz \
     $B/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz
curl -o data/geuvadis/EUR373.gene.cis.FDR5.best.rs137.txt.gz \
     $B/EUR373.gene.cis.FDR5.best.rs137.txt.gz
# sample overlap: GEUVADIS columns intersected with the phase-3 EUR panel
curl -s http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel \
  | awk '$3=="EUR"{print $1}' | sort > data/geuvadis/eur_phase3.txt
gzcat data/geuvadis/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz | head -1 \
  | tr '\t' '\n' | grep -E '^HG|^NA' | sort > data/geuvadis/geuvadis_samples.txt
comm -12 data/geuvadis/eur_phase3.txt data/geuvadis/geuvadis_samples.txt \
  > data/geuvadis/geuvadis_eur_overlap.txt   # 358 individuals
```

Headline result (best selected variant; **bold** tags the published lead at
r2 = 1.00; "---" = empty selection at the method's prespecified threshold):

| Gene | Lead eQTL | CBF-LMM | SuSiE | BayesB | BayesR | LMM scan |
|---|---|---|---|---|---|---|
| ERAP2 | rs2910686 | **rs2927608** | --- | **rs2910686** | --- | --- |
| RPS26 | rs10876864 | **rs10876864** | --- | **rs10876864** | **rs10876864** | --- |
| SLFN5 | rs11080327 | **rs11080327** | **rs11080327** | rs883416 (0.97) | **rs11080327** | **rs11080327** |
| SNHG5 | rs1059307 | **rs1059307** | **rs1059307** | **rs1059307** | **rs1059307** | --- |
| FLVCR1-AS1 | rs12123978 | **rs61832055** | --- | **rs11120042** | --- | --- |
| PEX6-region | rs6907751 | rs9986447 +1 (0.88) | **rs6907751** +1 | rs2296804 (0.63) | --- | --- |
| TRA2A-AS | rs10233039 | **rs6461691** | --- | rs10266123 (0.99) | --- | --- |
| ZNF266 | rs10420709 | **rs11878970** | --- | **rs10411141** | --- | --- |

CBF-LMM returns a selection at all eight loci and tags the published lead at
r2 = 1.00 at seven of eight; SuSiE's PIP > 0.99 selection is empty at five
loci (posterior mass split across perfect-LD proxies).


Each script supports `--cores N` for `mclapply` parallelism and writes per-cell
RDS checkpoints; re-running an interrupted job resumes from the last completed
checkpoint.  **Important:** set
`OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1
MKL_NUM_THREADS=1` to avoid Accelerate/OpenBLAS thread contention under
`mclapply` workers (the most common source of order-of-magnitude slowdowns).

> Note: the comparison is restricted to methods that model the polygenic
> as sparse whole-genome regressions (BayesB, BayesR) or via credible-set
> decomposition (SuSiE), plus a single-kernel frequentist LMM scan. 

## Internal validation

A six-check validation suite is provided in
[`sim/validate_lmm_bayes.R`](sim/validate_lmm_bayes.R).  V1 confirms
machine-precision recovery of Wakefield's ABF at $\delta \to 0$; V2 checks
asymptotic agreement with plug-in REML at $n = 1{,}000$; V3–V6 cover type-I
behaviour under the null, quadrature convergence, stepwise recall, and the
joint-vs-marginal distinction.

```bash
Rscript sim/validate_lmm_bayes.R
```

## Installation

See [`INSTALL.md`](INSTALL.md).  Briefly:

```r
install.packages(c("statmod", "Matrix",                 # core
                   "dplyr", "tidyr", "ggplot2",         # tables/figures
                   "susieR", "rrBLUP", "BGLR",          # comparators + mouse data
                   "hibayes"))                          # BayesR comparator
```

All packages are on CRAN; no GitHub-only dependencies.

## Citation

If you use this code in your work, please cite the manuscript (see
[`CITATION.cff`](CITATION.cff)):

```
@article{Doulabe2026LMMBayes,
  author  = {Doulabe, Kossi N. and Lakhal-Chaieb, Lajmi},
  title   = {Conditional Bayes-Factor Selection with an Adaptive
             Polygenic Background},
  journal = {(submitted)},
  year    = {2026}
}
```

## License

MIT — see [`LICENSE`](LICENSE).
