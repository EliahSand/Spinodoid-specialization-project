# Spinodoid-specialization-project

Workflow and tooling for generating spinodoid unit cells, tiling them into larger lattices, and exporting watertight STL meshes for downstream CAD/FEA studies.

## MATLAB generators

- `Matlab/main.m` – entry point for creating a single periodic spinodoid cell (cubes).
- `Matlab/tile_spinodoid.m` – tiles the cube along X/Y/Z to create multi-cell lattices.
- `Matlab/make_spinodoid_sheet.m` – dense base slab + closed spinodoid layer (panel studies).
- `Matlab/make_spinodoid_relief_sheet.m` – plate with a constant-thickness spinodoid “carpet”.
- `Matlab/make_spinodoid_stent.m` – wraps the field onto a cylindrical shell for individual stents.
- `Matlab/tile_spinodoid_stent.m` – fuses several stents along the axial direction (stent strings).

Each script writes its STL plus a `run_log.txt` into `Matlab/results/<subfolder>/run_<timestamp>` so every dataset stays reproducible (parameters, mesh stats, RNG seeds, etc.).

## Requirements

- MATLAB R2022b+ (parallel FFTs benefit from newer releases)
- Image Processing Toolbox (optional; only needed for connectivity pruning)  
- GPU support is automatic when a CUDA-capable device is present for `tile_spinodoid`.

## Results

Generated meshes land under `Matlab/results/` in generator-specific subfolders:

- `results/cells/` – single cubes (`main.m`)
- `results/tiled/` – cubic tilings (`tile_spinodoid.m`)
- `results/sheets/` – base + spinodoid sheet (`make_spinodoid_sheet.m`)
- `results/relief_sheets/` – plate + relief (`make_spinodoid_relief_sheet.m`)
- `results/stents/` – single stents (`make_spinodoid_stent.m`)
- `results/stents_tiled/` – axial stent chains (`tile_spinodoid_stent.m`)

Each run folder contains the STL and a matching log capturing parameters, mesh counts, and solid fractions.
