# GNN target pipeline — steps 5 to 8

## Overview

Steps 5–8 construct the supervised learning target from raw Abaqus midpoint output and
verify that the representation is valid before GNN training begins. The target is the
out-of-plane displacement field `u3(s)` sampled along the deformed midpoint line at 128
evenly-spaced arc-length positions `s ∈ [0, 1]`. This is optionally compressed to K=16
PCA coefficients. All steps communicate through `sample_id` as the sole join key. PCA is
fit exclusively on the training split so that val and test reconstruction errors are
genuine measures of generalization, not artefacts of basis construction.

---

## Step 5 — Build `u3(s)` target dataset

**Script:** `Matlab/GNN/pipeline/step5_build_u3_target_dataset.m`

**Inputs (per sample):** `midpoint_results_shell.csv` (columns Y, Z, U2, U3), `meta.json`

**Outputs:**
- Per sample: `target_u3_profile.mat`
- Aggregated: `data/dataset/targets/u3_targets.mat` — `U3_mat` (M×128), `U2_mat`,
  `s_grid`, `sample_ids`, `tr_ratio`, `ang_deg`
- Log: `logs/step5_u3_target_log.txt`

### Design choices and reasoning

**Arc length computed on the deformed curve `(Y+U2, Z+U3)`**
The target describes the physical shape of the midline under load. Parameterizing by
deformed arc length means the grid spacing is physically uniform in the loaded state,
which is the configuration the GNN needs to predict. Using undeformed arc length would
distort the function being learned (high-curvature zones would be under-resolved).

**Node ordering by undeformed Y**
Before computing arc length, nodes are sorted by their undeformed Y coordinate. This
gives a stable, monotone ordering that is independent of the deformation magnitude and
direction. If ordering were done on the deformed Y (i.e., Y+U2), large in-plane
displacements could produce a non-monotone ordering and break the arc-length integral.

**`pchip` interpolation onto a 128-point fixed grid**
`pchip` (piecewise cubic Hermite) is shape-preserving and monotone between data points.
Standard `spline` can overshoot near kinks, which in the spinodoid midpoint profiles
would introduce spurious oscillations. 128 grid points (`N_GRID=128`) is dense enough
to resolve short-wavelength wrinkles in lamellar sheets while keeping the raw target
vector a manageable size before PCA.

**Hard QC gates**
Five failure modes are checked and logged explicitly rather than silently dropped:

| Code | Condition |
|---|---|
| `fail:missing_csv` | no Abaqus output file |
| `fail:nan_inf` | any NaN or Inf in Y, Z, U2, U3 |
| `fail:too_few_points(n)` | fewer than 64 raw midplane nodes |
| `fail:u3_out_of_range` | `max(|u3|) > 0.5 m` |
| `fail:non_monotone_arc_length` | degenerate geometry post-deformation |

These gates match the sample acceptance criteria in CLAUDE.md. A silent skip would
contaminate `U3_mat` with unreliable rows that could then enter the train split.

**PCA is not fit in this step**
Step 5 is split-agnostic — it does not know which samples will be training data. Fitting
PCA here would force a later re-run whenever the split changes and, more importantly,
would make it impossible to prevent val/test information from entering the basis.

---

## Step 6 — Stratified train / val / test split

**Script:** `Matlab/GNN/pipeline/step6_split_dataset.m`

**Inputs:** `u3_targets.mat`

**Outputs:** `data/dataset/targets/split_indices.mat` — `train_idx`, `val_idx`,
`test_idx`, `rng_seed`, `train_frac`, `val_frac`, `test_frac`

### Design choices and reasoning

**Stratification on `(tr_ratio, ang_deg)` bins**
Samples in the dataset are generated in increasing lamellar order within each
`(tr_ratio, ang_deg)` parameter combination. A naive contiguous 70/15/15 split would
assign all `run01`–`run07` samples to train, leaving val and test without representation
of certain geometry classes. Stratification splits each bin independently so every
parameter combination appears in all three partitions.

**70 / 15 / 15 split fractions**
The 70 % training share maximizes the size of the basis fitted in step 7. The 15 % val
share is large enough for PCA reconstruction diagnostics (step 8) to be reliable across
hundreds of samples. Test is held out entirely until after hyperparameter decisions are
finalized.

**Minimum 1 / 1 / 1 per bin with test-priority fallback**
For bins with very few samples (< 5), the fallback prioritizes giving at least one
sample to test and one to val before assigning the remainder to train. This ensures no
split is completely empty for any parameter combination.

**Deterministic `rng(42)`**
Using a fixed seed means the split is reproducible and the seed is recorded in the
manifest. The Q3 stability check in step 8 verifies the PCA basis is insensitive to an
alternate seed (43), confirming reproducibility is not masking a fragile split.

---

## Step 7 — Fit PCA on training data only

**Script:** `Matlab/GNN/pipeline/step7_fit_pca_u3.m`

**Inputs:** `u3_targets.mat`, `split_indices.mat`

**Outputs:**
- `pca_model.mat` — `coeff` (128×16), `u3_mean` (1×128), `explained`, `n_components`
- `pca_targets.mat` — `Z_train`, `Z_val`, `Z_test` (each n×16)

### Design choices and reasoning

**PCA fit on training data only**
If the mean and basis were computed from all M samples, val and test profiles would
already be partially encoded in the subspace. The resulting reconstruction errors would
underestimate how well the basis generalizes to unseen data. Fitting on train only
makes the step 8 val/test RMSE a genuine out-of-sample measure.

**K = 8 components**
The scree plot (step 8 Q2) shows that 5 components already reach 99 % variance and 9
reach 99.9 %. K=8 is therefore chosen as a round number that sits above the 99 %
threshold with a margin of 3 components, while halving the GNN output size relative to
the original K=16. The GNN has 8 regression targets instead of 128, and the
reconstruction round-trip loses less than ~0.05 % of variance on average.

**Mean-centering with the training mean**
Both the mean and the basis are derived from train only. Val and test are projected as
`Z = (U3 − u3_mean_train) × coeff`. Using a per-split mean would break the consistent
affine transform and introduce an implicit data-normalization step that the GNN would
need to learn to invert.

**Note on `explained` from MATLAB's `pca`**
MATLAB's `pca(X, 'NumComponents', K)` returns `explained` as the variance spectrum for
*all* components (full rank of the training matrix), not just the K requested ones.
Downstream code that uses `explained` must therefore slice to `explained(1:K)` before
comparing against `coeff` — see the fix in step 8 line 104.

---

## How PCA encodes and reconstructs deformed profiles

### What PCA finds

Each sample produces a 128-point displacement vector `u3 ∈ ℝ¹²⁸`. Across the training
set this gives a matrix `U3_train` of shape `(n_train × 128)`. PCA finds K orthogonal
directions in this 128-dimensional space — the **principal components** — ordered so
that the first direction captures the most variance, the second the most of what remains,
and so on. These directions are the columns of `coeff` (128×K).

Geometrically, the principal components are the "shapes" that most commonly appear in
the dataset. PC1 for a sheet under bending is typically a smooth global bow — the
dominant deformation mode. PC2 might be an anti-symmetric twist. PC3 onwards capture
progressively finer local wrinkle patterns.

### Encoding: 128 numbers → K numbers

Given a new profile `u3` (1×128), the encoding step is:

```
z = (u3 − u3_mean) × coeff          % shape: 1 × K
```

`u3_mean` is the mean profile of the training set. Subtracting it centers the data so
the PCA axes describe *deviations* from the average shape. Each element `z_k` is the
**dot product** of the mean-subtracted profile with the k-th principal component — it
measures how strongly that mode appears in this particular deformation. A large positive
`z_1` means the sample bends strongly in the direction of the dominant training mode; a
near-zero `z_1` means it barely participates in that mode.

This projection is lossless only if the profile lies exactly in the K-dimensional
subspace. In practice, the discarded dimensions (components K+1 … 128) hold the
residual that cannot be represented.

### Decoding: K numbers → 128 numbers

The reconstruction is:

```
u3_recon = z × coeff' + u3_mean      % shape: 1 × 128
```

This is a **weighted sum of the K basis shapes**, plus the mean profile:

```
u3_recon(s) = u3_mean(s)  +  z_1 · coeff(:,1)  +  z_2 · coeff(:,2)  + … +  z_K · coeff(:,K)
```

Each basis shape `coeff(:,k)` is a fixed 128-point curve (a spatial mode). The
coefficient `z_k` scales how much of that mode is added back. The reconstruction is
therefore a linear combination of K pre-learned shapes — the GNN's job is to predict
the K scalar weights `z` given the structural graph as input.

### Why reconstruction is lossy and what the error means

The truncation from 128 to K dimensions discards the components that carry the least
variance *on average*. For most samples this residual is negligible. But some samples
have displacement profiles with fine spatial detail (short-wavelength wrinkles, sharp
kinks) that are not well-captured by the dominant modes. For those samples the
reconstruction error can be significantly higher than the dataset average.

This is why the 99 % average-variance threshold does not guarantee low worst-case error.
A sample whose deformation is mostly in component 6 will reconstruct poorly with K=5
even though that component is "small" across the training set as a whole.

**What the GNN actually learns:** the GNN predicts `z` (K numbers) from the graph. At
inference, `z` is decoded back to `u3_recon` using the fixed saved basis. The
reconstruction error sets a lower bound on GNN accuracy — even a perfect GNN cannot
predict the discarded residual. Choosing K is therefore a trade-off:

| K too small | K too large |
|---|---|
| High reconstruction ceiling — even perfect coefficients give bad `u3` | More GNN output targets — harder regression, more data needed |
| Worst-case samples dominated by residual | Targets include near-noise components |

A practical rule: choose the smallest K such that the worst-case reconstruction error
(from step 8 Q5 **val** plots) is below the accuracy you need from the full pipeline.

### Why worst cases are shown on the validation set

The basis is fit on training data, so training reconstruction is always optimistic —
the basis was built to describe those exact samples. Validation samples were held out
during PCA fitting, so their reconstruction error is an honest measure of how well
the K-dimensional subspace generalizes to unseen geometry. The worst val samples are
the shapes farthest from the training subspace — showing them is a gut-check for
whether those samples should be in the dataset at all, or whether K needs increasing.

---

## Step 8 — Diagnostics (Q1–Q6)

**Script:** `Matlab/GNN/pipeline/step8_verify_pca_target.m`

**Inputs:** `u3_targets.mat`, `split_indices.mat`, `pca_model.mat`, `pca_targets.mat`

**Outputs:**
- `diagnostics/pca_diagnostics.json`
- `diagnostics/scree.png`
- `diagnostics/coeff_hist.png`
- `diagnostics/worst_recon/worst_NN_<sample_id>.png` (one file per worst val sample)

### The six checks and their purpose

**Q1 — Reconstruction fidelity**
PCA is a lossy compression: the K-dimensional coefficients Z are computed as
`Z = (U3 − u3_mean) × coeff`, and the profile is recovered as
`U3_recon = Z × coeff' + u3_mean`. Reconstruction error measures how much of the
original 128-point shape is lost in that round-trip. Step 8 computes this for every
sample in each split and reports:

- *Absolute RMSE* — mean over samples of `sqrt(mean((u3_true − u3_recon)²))` in metres.
  Gives a physical sense of the error scale.
- *Relative RMSE* — absolute RMSE divided by `max(|u3_true|)` per sample, then averaged.
  Normalizes for the fact that some samples deflect much more than others.

The pass condition is `val_rel / train_rel < 2.0` and `test_rel / train_rel < 2.0`.
A ratio near 1.0 means the basis generalizes well: the components learned from the
training profiles describe val/test shapes equally well. A large ratio (> 2) would
indicate the training set does not span the full shape space, and the PCA basis is
essentially extrapolating for some val/test samples. The train relative RMSE being `Inf`
in the first run was caused by samples with near-zero displacement — those have
`max(|u3|) ≈ 0`, making the denominator unstable; such samples should be flagged by the
`u3_out_of_range` gate in step 5 or handled with a minimum-displacement filter.

**Q2 — Scree / variance compactness**
Reports the number of components required to reach 90 / 95 / 99 / 99.9 % cumulative
variance. Passes if K=16 is sufficient for 99 %. This is a fast sanity check on the
choice of K — if the profile were unexpectedly complex (e.g., after changing geometry
parameters), the check would flag that more components are needed.

**Q3 — Basis stability**
Refits PCA on an alternate-seed (43) stratified train split and measures the subspace
angle between the two K-dim bases plus the mean per-component cosine similarity. Passes
if the angle is < 10°. A large subspace angle would mean the basis is a statistical
accident of the particular seed-42 sample ordering, which would make the GNN coefficients
incomparable across experiments.

**Q4 — Coefficient distributions**
Histograms of Z scores (train) for each PC with std and skewness. Not a hard pass/fail.
Flags non-Gaussian or bimodal distributions that could make gradient-based optimization
unstable. PC1 typically carries the global bending mode and should be approximately
Gaussian given the range of geometries sampled.

**Q5 — Worst reconstructions**
Identifies the top-5 val samples by relative reconstruction error and saves each as an
individual PNG (`worst_recon/worst_01_<sample_id>.png` etc.). Individual files are used
instead of a combined panel so that axis labels, titles, and profile shapes are
readable at native resolution. These are a visual gut-check — if the worst cases look
like reasonable shapes that merely have high-frequency content, the PCA is working as
intended. If they look like failed simulations, the QC gates in step 5 need tightening.

**Q6 — Learnability baseline**
Fits ridge regression from four coarse scalar inputs (`tr_ratio`, `ang_deg`,
`sin(ang_deg)`, `cos(ang_deg)`) onto each of the first four PC scores and reports R² on
val. If R² > 0 the bulk geometry parameters explain some variance — expected. If R² < 0
the PC score is noisier than the mean prediction, which would be a red flag. The intent
is the opposite of a strong baseline: we *want* the GNN to need graph-level structure to
do better than these four scalars.

**Why all checks are warnings, not hard errors**
The thresholds (relative RMSE ratio < 2, subspace angle < 10°, k_99 ≤ 32) are
informed guesses. Failing one check does not necessarily mean the dataset is unusable —
it means the author should read the `pca_diagnostics.json` and the plots before deciding
whether to retrain or adjust K. Hard errors here would block valid datasets on
borderline cases.

---

## How the steps fit together

```
data/dataset/samples/*/midpoint_results_shell.csv
          |
          v step 5
data/dataset/targets/u3_targets.mat   (M × 128 matrix, all QC-passed samples)
          |
          v step 6
data/dataset/targets/split_indices.mat  (train/val/test index vectors)
          |
          v step 7  (reads both files above; fits on train rows only)
pca_model.mat   (128×16 basis + train mean)
pca_targets.mat (Z_train, Z_val, Z_test — n×16 coefficient matrices)
          |
          v step 8
diagnostics/pca_diagnostics.json + *.png
```

---

## Join key discipline

Every step indexes into aggregates by `sample_id`. `u3_targets.mat` preserves the row
order established in step 5, and `split_indices.mat` contains integer indices into that
same row ordering. `pca_targets.mat` stores `train_ids`, `val_ids`, `test_ids` alongside
the coefficient matrices so that any downstream script can verify alignment without
relying on row position alone. Never match samples by filename pattern when `sample_id`
exists.
