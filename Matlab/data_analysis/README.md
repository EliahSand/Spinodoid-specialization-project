# Spinodal solid–shell analysis

- Entry point: `data_prep/run_all_analysis.m`
- Scope is intentionally minimal:
  - overlay plots (solid vs shell)
  - curvature comparisons (raw overlay + baseline-offset removed + theory check)
- Generated outputs are written under `results/analysis/single_runs/`.
- Run discovery is automatic under `../results/sheets/**` for folders that contain both `FEA` and `FEA_shell`.
