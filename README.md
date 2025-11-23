# Spinodoid-specialization-project
Author: Eliah Sand (NTNU, 2025)

End-to-end pipeline for generating periodic spinodoid lattices (cells, sheets, stents), exporting meshes, and running Abaqus solid/shell simulations with mid-plane probes.

## Theory (periodic GRF)

The spinodal field is synthesized in Fourier space as a periodic Gaussian random field:
$$
\phi(\mathbf{x}) = \sum_{m=1}^{M} a_m \cos\!\left( \frac{2\pi}{L}\, \mathbf{k}_m \!\cdot\! \mathbf{x} + \gamma_m \right),
$$
with $\mathbf{k}_m = (i,j,k) \frac{2\pi}{L}$ and integer $(i,j,k)$. Every mode divides the box length $L$, so values and derivatives match on opposite faces:
$$
\phi(0,y,z) = \phi(L,y,z), \quad \partial_x \phi(0,y,z)=\partial_x \phi(L,y,z),
$$
and similarly for $y,z$. Stent variants build the field in $(\theta,z,r)$ with $k_\theta=n$ to enforce angular periodicity; padding/wrapping yields watertight cylindrical masks.

A percentile threshold converts $\phi$ to a binary mask; an isosurface of that mask becomes the STL. Periodicity → tileable geometry and clean PBCs.

## Generation workflow (MATLAB)

- Generators live under `Matlab/` (see Matlab/README.md for details and commands).
- Each run writes `results/.../run_<timestamp>/` containing:
  - STL
  - `run_log.txt` (parameters, thicknesses, spacing, solid fractions)
  - Mask `.mat` (logical volume + spacing/origin; sheets include `zVoxelThickness`)
  - `mesh_manifest.json` pointing to the mask/variable

## Mesh conversion (Python)

- Solids (voxel C3D8):  
  - Single: `python exp/mat_to_inp.py --mat <mask>.mat --var <var> --spacing <L/N>`  
  - Batch: `python exp/batch_mat_to_inp.py --roots Matlab/results --overwrite`
- Shells (S4R):  
  - Single: `python exp/mat_to_shell_inp.py --mat Matlab/results/.../sheet.mat --var sheetMask`  
  - Batch: `python exp/batch_mat_to_shell_inp.py --roots Matlab/results --overwrite`  
    Uses `t_base_mm` / `t_spin_mm` from `run_log.txt` when present.

## Abaqus runners (CAE noGUI)

- Solids: `abaqus cae noGUI=exp/run_spinodal_static.py -- <solid.inp>`  
  Batch: `python exp/batch_run_spinodal_solid.py --roots Matlab/results --overwrite`
- Shells: `abaqus cae noGUI=exp/run_spinodal_shell_static.py -- <shell.inp>`  
  Batch: `python exp/batch_run_spinodal_shell.py --roots Matlab/results --overwrite`

Outputs: ODB + mid-plane CSV/RPT under each run’s `FEA/` or `FEA_shell/`. Use `--dry-run` on batch tools to preview commands.

## Folder map

- `Matlab/` – generators, helpers, README (usage and file management)
- `exp/` – converters (`mat_to_inp.py`, `mat_to_shell_inp.py`), runners, batch tools
- `Matlab/results/` – generated datasets (STL, run_log.txt, mask .mat, manifest, FEA outputs)
