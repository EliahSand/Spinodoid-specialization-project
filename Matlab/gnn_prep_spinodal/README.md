# gnn_prep_spinodal

MATLAB module for converting Abaqus shell mesh + nodal CSV outputs into:

- a full reference graph
- a 1-node-thick skeleton graph
- a reduced structural graph

This module is now centered on the structural-graph workflow. Legacy generic downsampling methods are not part of the active pipeline.

## What It Does

- Parses Abaqus `.inp` files:
  - node labels and coordinates (`*Node`)
  - element connectivity (`*Element`)
  - element sets (`*Elset`, including `generate`)
- Selects a region from an elset (for example `SPINODAL_SHELL`) or auto-detects a spinodal/top elset.
- Builds a full reference graph:
  - nodes = selected mesh nodes
  - edges = unique shell edges derived from element connectivity
  - boundary nodes from edge use-count (edges used by only one element)
- Maps nodal CSV fields (`Label`, `U1`, etc.) to graph node indices.
- Extracts a 1-node-thick skeleton/centerline graph from the full graph.
- Compresses the skeleton into a reduced structural graph while preserving topology and path geometry.
- Exports graph tables and plots.

## Folder Layout

```text
gnn_prep_spinodal/
  src/
    graph/
    io/
    mesh/
    vis/
  examples/
  tests/
  data/
  out/
```

## Quick Start

Run the hardcoded entry script:

```matlab
run_main_spinodal_gnn_prep
```

Or call the main function directly:

```matlab
outputs = main_spinodal_gnn_prep(inpPath, csvPath, ...
    'ElsetName', 'SPINODAL_SHELL', ...
    'StructuralDetailLevel', 0.20, ...
    'OutDir', 'Matlab/gnn_prep_spinodal/out');
```

Minimal example script:

```matlab
run('Matlab/gnn_prep_spinodal/examples/example_build_graphs.m')
```

## Main Output

- `full_graph_*`: full reference graph export
- `skeleton_graph_*`: extracted 1-node-thick skeleton graph export
- `structural_graph_*`: reduced structural graph export
- `*.mat`: serialized graph structs
- optional `.png` plots for:
  - full reference graph
  - structural overlay

## Out Folder Structure

Each run is written to its own subfolder:

```text
out/
  <prefix>_<timestamp>/
    metadata/run_info.txt
    full_graph/
      data/
      plots/
    structural_graph/
      data/
      plots/
```

## Structural Graph

Build the reduced skeleton/centerline graph directly from the dense reference graph:

```matlab
outputs = main_spinodal_gnn_prep(inpPath, csvPath, ...
    'ElsetName', 'SPINODAL_SHELL', ...
    'StructuralDetailLevel', 0.20, ...
    'StructuralMinIslandNodes', 1);
```

This produces:

- `skeleton_graph_*`: 1-node-thick extracted skeleton graph
- `structural_graph_*`: compressed structural graph with:
  - reduced node coordinates
  - reduced edges
  - edge path length
  - edge polyline
  - edge-to-full-node-path mapping
  - node mapping back to dense graph node ids
- `structural_graph_overlay.png`: dense + skeleton + reduced graph overlay

## Graph Layers

The structural overlay plot shows several representations of the same spinodal morphology.

### Dense

This is the original reference graph built directly from the selected Abaqus shell region.

- nodes = all selected mesh nodes
- edges = local mesh connectivity from the shell elements

This is the highest-resolution graph and acts as the source geometry/topology for all later steps.

### Skeleton

This is the extracted 1-node-thick centerline network.

- the dense XY occupancy is rasterized on the native node grid
- the occupied region is thinned with a topology-preserving skeletonization step
- the result is converted back into a graph

This graph captures the medial paths, loops, branches, endpoints, and junctions of the morphology.

### Skeleton Nodes

These are the node locations of the skeleton graph.

- they are the discrete centerline sample points after thinning
- they often still include many degree-2 chain nodes along long branches

They are shown separately so you can see how finely the skeleton is sampled.

### Structural

This is the reduced graph built from the skeleton graph.

- degree-2 chains are collapsed into single reduced edges
- junctions and endpoints are preserved
- extra waypoints can be retained depending on `StructuralDetailLevel`

This is the simplified graph intended to represent the main structural morphology with much lower complexity.

Each reduced edge still carries the original path information through:

- `edge_length`
- `edge_polyline`
- `edge_full_node_paths`

So a structural edge is a compressed representation of a real path along the skeleton, not just a straight shortcut.

### Structural Nodes

These are the reduced graph nodes kept after chain compression.

They correspond to important anchor points such as:

- endpoints
- junctions
- optional intermediate waypoints when more detail is requested

Lower `StructuralDetailLevel` gives fewer structural nodes and a coarser graph. Higher `StructuralDetailLevel` keeps more waypoints and produces a finer reduced graph.

## Plot Tuning

Useful options for clearer plots:

```matlab
'PlotNodeMarkerSize', 3, ...
'PlotBoundaryMarkerSize', 3.5, ...
'PlotBoundaryMode', 'interior', ...
'PlotUnitsScale', 1, ...
'PlotXLabel', 'x [mm]', ...
'PlotYLabel', 'y [mm]'
```
