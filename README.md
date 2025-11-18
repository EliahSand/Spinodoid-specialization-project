# Spinodoid-specialization-project

Workflow and tooling for generating spinodoid unit cells, tiling them into larger lattices, and exporting watertight STL meshes for downstream CAD/FEA studies.

## MATLAB generators

- `Matlab/main.m` – entry point for creating a single periodic spinodoid cell (cubes).
- `Matlab/main_cubed.m` – formerly the tiler; replicates the `main` cube along X/Y/Z and fuses seams.
- `Matlab/PSSH.m` – *Periodic Spinodal Sheet Layered*: dense base slab with an extruded spinodal layer.
- `Matlab/PSS.m` – *Periodic Spinodal Stent*: cylindrical spinodal stent with angular periodization (`theta_partitions`) and optional axial tiling (`tilesAxial`).
- `Matlab/PSSL.m` – *Periodic Spinodal Stent Layered*: dense inner wall + spinodal relief outer layer, sharing the same periodization controls.

Each script writes its STL plus a `run_log.txt` into `Matlab/results/<subfolder>/run_<timestamp>` so every dataset stays reproducible (parameters, mesh stats, RNG seeds, etc.).

## Requirements

- MATLAB R2022b+ (parallel FFTs benefit from newer releases)
- Image Processing Toolbox (optional; only needed for connectivity pruning)  
- GPU support is automatic when a CUDA-capable device is present for `tile_spinodoid`.

## Results

Generated meshes land under `Matlab/results/` in generator-specific subfolders:

- `results/cells/` – single cubes (`main.m`)
- `results/tiled/` – cubic tilings (`main_cubed.m`)
- `results/sheets/` – layered sheets (`PSSH.m`)
- `results/stents/single|periodic/<type>/` – spinodal stents classified by whether `theta_partitions` or `tilesAxial` exceed 1 (`PSS.m`)
- `results/stent_relief/single|periodic/<type>/` – layered stents, likewise split into single vs periodized (`PSSL.m`)

Each run folder contains the STL and a matching log capturing parameters, mesh counts, and solid fractions.

### Voxel meshing for Abaqus

If you need a voxelized C3D8 mesh for Abaqus without installing extra Python packages, export the MATLAB mask (e.g., `solidMask`, `stentMask`) to a `.mat` file and run:

```bash
python exp/mat_to_inp.py --mat solidMask.mat --var solidMask --spacing <L/N> \
    --peel 3 --output exp/results/spinodal_mesh.inp --density 1100 --elastic 2.5e9 0.3
```

The script depends only on NumPy + SciPy (plus h5py for v7.3 .mat files) and produces an Abaqus-ready voxel mesh.
