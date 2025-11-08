# Spinodal Periodic Structure Generator

This folder contains several MATLAB generators tailored for periodic spinodoids:

- `main.m` (single cell) – produces one periodic spinodal microstructure suitable for additive manufacturing, finite element homogenization, and metamaterial design.  
- `tile_spinodoid.m` (multi-cell tiler) – reuses the single-cell recipe, tiles it along X/Y/Z, and exports a much larger seamless lattice with its own run log.  
- `make_spinodoid_sheet.m` – stacks a dense base slab with a closed spinodoid layer cut from the same periodic field (good for panel studies).  
- `make_spinodoid_relief_sheet.m` – fuses a homogeneous plate with a constant-thickness spinodoid “carpet” extruded from a 2D slice, matching the periodic-with-plate style reference.
- `make_spinodoid_stent.m` – wraps the periodic spinodal field onto a cylindrical shell (Ri/Ro/H) to export a watertight single stent using the same padding/wrapping workflow as the reference stent generator.
- `tile_spinodoid_stent.m` – replicates the stent volume along the axial direction and fuses the segments into one continuous lattice, analogous to how `tile_spinodoid.m` extends the cubic cells.

All shared utilities live under `helpers/` and are added to the MATLAB path by each script at runtime.

Stent scripts (`make_spinodoid_stent.m` and `tile_spinodoid_stent.m`) rely on cylindrical padding helpers embedded in the files, so you can run them without any external toolboxes: just `cd Matlab`, call the script, and check `results/stents/` or `results/stents_tiled/` for the STL/log pair.

Both scripts emit watertight STL solids with controlled porosity, anisotropy, and feature size.

---

## 1. Overview

The generator builds a 3D scalar field $ \phi(\mathbf{x}) $ representing a Gaussian random field (GRF) synthesized directly
in Fourier space.  
A percentile threshold converts $ \phi(\mathbf{x}) $ into a binary mask defining the solid region.  
An isosurface of this binary mask is exported as a watertight STL file.

Unlike conventional spinodal-like generators, this method ensures **true periodicity** by mathematical construction —  
both function values and derivatives match perfectly at opposite faces of the cube.  
This makes the structure tileable and suitable for simulations using periodic boundary conditions (PBC).

---

## 2. What makes it periodic?

The periodicity comes from how the scalar field $ \phi(\mathbf{x}) $ is defined:

$$
\phi(\mathbf{x}) \;=\; \sum_{m=1}^{M} a_m \,\cos\!\left( \frac{2\pi}{L}\, \mathbf{k}_m \!\cdot\! \mathbf{x} \;+\; \gamma_m \right),
$$

where

- $ L $ — box length  
- $ \mathbf{k}_m = (i, j, k)\, \dfrac{2\pi}{L} $, with $ i, j, k \in \mathbb{Z} $ (integers)  
- $ \gamma_m $ — random phase for each mode

Because every cosine term has a wavelength that divides the box length $ L $ exactly,
the sum $ \phi(\mathbf{x}) $ is **analytically periodic**:

$$
\begin{aligned}
\phi(0, y, z) &= \phi(L, y, z),\\[4pt]
\frac{\partial \phi}{\partial x}(0, y, z) &= \frac{\partial \phi}{\partial x}(L, y, z),
\end{aligned}
$$

and the same holds for $ y $ and $ z $.  
Thus, opposite faces of the cube wrap seamlessly — both in value and in slope.

When you threshold $ \phi(\mathbf{x}) $ into a binary mask (solid vs. void), the interface defined by $ \phi = \text{const} $
also wraps continuously across boundaries, forming a **toroidal (periodic)** field.

---

## 3. Why this matters

- **Tiling** → Multiple unit cells can be placed side by side without visible seams  
- **Simulation** → Enables perfect PBCs in FEA or FFT homogenization  
- **Consistency** → No artificial edge effects  
- **Control** → You can analytically tune anisotropy, porosity, and feature wavelength

---

## 4. Parameters

| Parameter | Meaning | Typical values |
|------------|----------|----------------|
| `N` | Grid size (NxNxN). Higher = smoother geometry | 64–256 |
| `L` | Physical box length | 1.0 |
| `lambda_vox` | Target feature wavelength (voxels) | 16–32 |
| `bandwidth` | Spectral shell thickness (smaller = smoother) | 0.10–0.25 |
| `nModes` | Number of Fourier modes | 1000–10000 |
| `solid_frac` | Solid volume fraction (0–1) | 0.3–0.5 |
| `coneDeg` | Cone half-angles [θx θy θz] (90° = isotropic) | [90 90 90] |

---

## 5. Output

Each run writes a watertight STL file of a periodic spinodal solid:

```
spinodoid_N128_sf50_vox.stl
```

If a file with the same name already exists, the script automatically appends a version number:

```
spinodoid_N128_sf50_vox_v1.stl
spinodoid_N128_sf50_vox_v2.stl
```

### Tiled generator (`tile_spinodoid.m`)

The tiler follows the same periodic pipeline, then replicates the solid mask along each axis (`tiles = [tx ty tz]` or scalar). Key controls:

- `sigma_vox` – optional Gaussian blur applied in periodic Fourier space before thresholding.  
- `ups` – periodic upsampling factor (1 = off, 2/3 = smoother STL).  
- `tiles` – replication counts per axis; accepts scalar shorthand.  
- `skin_thickness_vox` / `keep_largest_solid` – same closure/pruning controls as the single-cell script.  
- Automatic CPU multi-threaded FFTs, with opportunistic GPU execution when `gpuDeviceCount > 0`.

All scripts now drop artefacts into `Matlab/results/` but separated by generator for easy navigation:

- `results/cells/` (`main.m`)
- `results/tiled/` (`tile_spinodoid.m`)
- `results/sheets/` (`make_spinodoid_sheet.m`)
- `results/relief_sheets/` (`make_spinodoid_relief_sheet.m`)
- `results/stents/` (`make_spinodoid_stent.m`)
- `results/stents_tiled/` (`tile_spinodoid_stent.m`)

Inside each subfolder, every run gets its own timestamped directory (e.g., `isotropic_N128_sf50/run_<timestamp>`), containing the STL and a concise `run_log.txt` that captures parameters, mesh statistics, and any post-processing steps. This keeps artefacts traceable even when many studies run in parallel.

---

## 6. Comparison with other generation methods

| Method | Periodicity | Computation | Notes |
|--------|--------------|-------------|-------|
| **This generator** | Exact (Fourier-based) | Instant | Tileable, analytical control |
| **Cahn–Hilliard PDE** | Approximate (numerical BCs) | Heavy | Physically realistic coarsening |
| **Random GRF (non-periodic)** | None | Instant | Edge discontinuities, non-tileable |

---

## 7. Smoothness recommendations

To increase smoothness:

1. **Increase resolution:** `N = 128` → `256`  
2. **Reduce spectral width:** `bandwidth = 0.12`  
3. **Apply Gaussian blur:** $\sigma \approx 1.0$ voxel (via FFT)  
4. **Upsample $\phi(\mathbf{x})$:** $2\times$ using periodic FFT zero-padding

All these steps operate in periodic space, preserving seamless boundaries.

---

## 8. Notes

- Exported surfaces are watertight and manifold.
- Filtering and smoothing are performed via FFT (periodic domain).
- The `skin_thickness_vox` parameter can add an external solid shell if desired.

---


### Author

Generated from **Eliah Sand’s specialization project (NTNU, 2025)**.  
If used academically, cite this as the *Periodic Spinodal Fourier Synthesis Method*.
