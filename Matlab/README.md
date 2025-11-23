# MATLAB Spinodal Generators

This folder hosts the spinodoid generators and their outputs; theory lives in the top-level README.

## Generators
- `main.m` – single periodic cell (cube)
- `main_cubed.m` – tiled cubes (X/Y/Z replication)
- `PSSH.m` – layered sheet (dense base + spinodal layer, periodic in XY)
- `make_spinodoid_relief_sheet.m` – relief sheet variant
- `make_spinodoid_stent.m` / `tile_spinodoid_stent.m` – cylindrical stents (single or tiled axially)

Helpers sit in `helpers/` and are added to the MATLAB path by each script.

## Outputs per run
Each run writes a timestamped folder under `results/.../` containing:
- STL
- `run_log.txt` (parameters, thicknesses, spacing, solid fractions)
- Mask `.mat` (logical array + spacing/origin; sheets include `zVoxelThickness`)
- `mesh_manifest.json` (points to the mask/variable)

## Typical usage
```matlab
cd Matlab
PSSH;              % layered sheet
make_spinodoid_stent;  % cylindrical stent
```
Inspect the generated `results/.../run_*` for STL/log/mask.

## Hand-off to Python/Abaqus
- Convert masks to meshes:
  - Solids: `python exp/mat_to_inp.py ...` or `python exp/batch_mat_to_inp.py --roots Matlab/results`
  - Shells: `python exp/mat_to_shell_inp.py ...` or `python exp/batch_mat_to_shell_inp.py --roots Matlab/results`  
    (reads `t_base_mm` / `t_spin_mm` from `run_log.txt` when present)
- Run FEAs:
  - Solids: `abaqus cae noGUI=exp/run_spinodal_static.py -- <solid.inp>`
  - Shells: `abaqus cae noGUI=exp/run_spinodal_shell_static.py -- <shell.inp>`
  - Batch runners: `exp/batch_run_spinodal_solid.py`, `exp/batch_run_spinodal_shell.py`
- Results land under each run’s `FEA/` or `FEA_shell/` (ODB + mid-plane CSV/RPT).
