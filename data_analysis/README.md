# SOLID vs SHELL comparison

Python tools for comparing SOLID and SHELL FEA results per node, focused on displacements (U1–U3), stress components (Sij), and von Mises (SMises). The scripts generate publication-quality plots and metrics with a predictable output layout.

## Expected data layout
- Default data root: `Matlab/results/sheets/lamellar`
- Each case should contain `FEA_shell/midplane_results_shell.csv`
- Matching SOLID results are expected in `FEA_solid/` next to the shell folder (e.g. `FEA_solid/midplane_results_solid.csv`). If SOLID data is missing, that case is skipped.
- Required columns in both CSVs: `Label`, `U1`, `U2`, `U3`, `S11`, `S22`, `S33`, `S12`, `S13`, `S23`, `SMises`.

## Output structure
All outputs go to `data_analysis/results` by default (override with `--output-root`):
- `results/metrics/<group>_metrics.csv` — MAE, RMSE, relative errors, and Pearson r per component and case.
- `results/plots/<group>/<case>/...` — per-case scatter (SOLID vs SHELL), relative error histograms, and absolute error boxplots.
- `results/plots/<group>/aggregate/...` — heatmaps and barplots summarizing errors/correlations across cases.

## Dependencies
- Python 3.10+
- `pandas`, `numpy`, `matplotlib`, `seaborn`

Install (if needed):
```bash
pip install pandas numpy matplotlib seaborn
```

## How to run
Run with `python -m` so package imports resolve:
```bash
python -m data_analysis.dataPrep.compare_displacements [--data-root PATH] [--output-root PATH]
python -m data_analysis.dataPrep.compare_stresses [--data-root PATH] [--output-root PATH]
python -m data_analysis.dataPrep.compare_von_mises [--data-root PATH] [--output-root PATH]
```
Useful flags:
- `--limit N` — process only the first N cases (quick smoke tests).

## What the scripts do
For each case with SOLID+SHELL CSVs:
- Merge results on `Label` and compute absolute and relative (percent) errors.
- Plot SOLID vs SHELL scatter with y=x reference and correlation `r`.
- Plot relative error distributions and absolute error boxplots.
- Export metrics (MAE, RMSE, max error, mean |relative error|, Pearson r).
- Build aggregate heatmaps/barplots across cases for quick comparison.
