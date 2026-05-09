# GNN Dataset Pipeline (Simple)

This pipeline reuses the existing scripts that already work.

No extra manifest system is used.  
Data flow is step-by-step with files stored directly under `Matlab/GNN/data`.

## Data Layout

```text
Matlab/GNN/data/
  raw/
    samples/      % PSSCone run folders + Abaqus outputs (in-place)
  dataset/
    samples/      % one folder per run_name with structural graph (.mat) + midpoint CSV
    splits/
```

## Steps

1. Generate 6000 shell geometries (91 angles, 3 ratios, unique seeds per angle-ratio):
   - `matlab -batch "run('Matlab/GNN/datasetCreation/step1_generate_psscone_shell_dataset.m')"`
2. Convert MAT -> shell INP:
   - `py -3 Matlab/GNN/datasetCreation/step2_mat_to_shell_inp.py`
3. Build structural graphs from INP only (pre-deformation):
   - `matlab -batch "run('Matlab/GNN/datasetCreation/step3_batch_structural_graph_from_inp.m')"`
4. Run Abaqus shell simulations and package midpoint CSVs:
   - `abaqus cae noGUI=Matlab/GNN/datasetCreation/step4_run_abaqus_shell.py --`
   - Dry-run: `abaqus cae noGUI=Matlab/GNN/datasetCreation/step4_run_abaqus_shell.py -- --dry-run`
   - Package only from existing `midplane_results_shell.csv` files (no Abaqus rerun):
     `abaqus cae noGUI=Matlab/GNN/datasetCreation/step4_run_abaqus_shell.py -- --package-only`

## Important policy

- Inside each `(angle, ratio)` group, `rngSeed` is never reused.
- `step3_batch_structural_graph_from_inp.m` writes one MATLAB file per run:
  `data/dataset/samples/<run_name>/sample.mat` (contains structural graph).
- Step 4 moves `FEA_shell/midplane_results_shell.csv` into:
  `data/dataset/samples/<run_name>/midpoint_results_shell.csv`.
- Step 4 does not duplicate graph CSV files.
- Abaqus output kept for targets is `midplane_results_shell.csv` (plus lightweight run metadata outside `FEA_shell`).
- `step4_run_abaqus_shell.py` does immediate cleanup of heavy Abaqus artifacts (including `*_job.inp` and `.jnl`) after successful CSV extraction.
