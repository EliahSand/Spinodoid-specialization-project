# gnn_prep_spinodal

MATLAB module for converting Abaqus shell mesh + nodal CSV outputs into graph data for GNN preprocessing.

## What It Does

- Parses Abaqus `.inp` files:
  - node labels and coordinates (`*Node`)
  - element connectivity (`*Element`)
  - element sets (`*Elset`, including `generate`)
- Selects a region from an elset (e.g., `SPINODAL_SHELL`) or auto-detects a spinodal/top elset.
- Builds a full reference graph:
  - nodes = selected mesh nodes
  - edges = unique shell edges derived from element connectivity
  - boundary nodes from edge use-count (edges used by only one element)
- Maps nodal CSV fields (`Label`, `U1`, etc.) to graph node indices.
- Builds a deterministic downsampled graph with tunable detail:
  - boundary retention controlled by `BoundaryKeepRatio`
  - interior sampling methods: `farthest`, `grid`, `hybrid`
  - optional shortest-path augmentation to preserve connectivity/topology
- Exports graph tables and simple plots.

## Folder Layout

```
gnn_prep_spinodal/
  src/
    io/
    mesh/
    graph/
    downsample/
    vis/
  examples/
  tests/
  data/
  out/
```

## Quick Start

Run:

```matlab
run('Matlab/gnn_prep_spinodal/examples/example_build_graphs.m')
```

Method comparison example:

```matlab
run('Matlab/gnn_prep_spinodal/examples/example_compare_downsample_methods.m')
```

The comparison script also runs an experimental midpoint graph generator from:
`src/downsample/downsample_midpoint_experimental.m`.

Or call directly:

```matlab
outputs = main_spinodal_gnn_prep(inpPath, csvPath, ...
    'ElsetName', 'SPINODAL_SHELL', ...
    'DetailLevel', 0.25, ...
    'DownsampleMethod', 'hybrid', ...
    'OutDir', 'Matlab/gnn_prep_spinodal/out');
```

## Main Output

- `*_nodes.csv`: node indices, labels, coordinates, boundary mask, optional node attributes.
- `*_edges.csv`: graph edge list with node indices + node labels.
- `*.mat`: serialized graph struct.
- Optional `.png` plots for full and downsampled graphs.

## Out Folder Structure

Each run is written to its own subfolder:

```
out/
  <prefix>_<timestamp>/
    metadata/run_info.txt
    full_graph/
      data/
      plots/
    downsampled_graph/
      data/
      plots/
    line_graph/            % only when DownsampleOutputMode='line'
      data/
      plots/
```

## Extra Downsampling (Single Line Graph)

Use:

```matlab
outputs = main_spinodal_gnn_prep(inpPath, csvPath, ...
    'DetailLevel', 0.08, ...
    'DownsampleOutputMode', 'line', ...
    'LineNodeBudget', 40);
```

This creates a component-wise line-like graph from the downsampled mesh.
`DownsampleOutputMode='line'` simplifies boundary components separately, so multiple major spinodal parts are represented (not just one global line).

To also downsample straight boundary segments in mesh mode:

```matlab
'BoundaryKeepRatio', 0.25, ...
'BoundaryStraightTolDeg', 6
```

When `BoundaryKeepRatio < 1`, boundary polylines are re-linked with shortcut edges
between kept boundary nodes so each spinodal boundary remains connected as a line.

## Compare Methods Quickly

Try different deterministic interior selection methods:

```matlab
'DownsampleMethod', 'farthest'   % geodesic spread from boundary seeds
'DownsampleMethod', 'grid'       % spatial bin representatives
'DownsampleMethod', 'hybrid'     % junction anchors + grid/farthest fill
```

For `hybrid`, tune:

```matlab
'HybridAnchorFraction', 0.35
```

## Experimental Midpoint Variant

Experimental-only function:

```matlab
midGraph = downsample_midpoint_experimental(downsampledMeshGraph, ...
    'UseInteriorBoundaryOnly', true, ...
    'ColumnStride', 5, ...
    'PairMode', 'extremes', ...
    'ComponentKeepRatio', 0.04, ...
    'MaxNodesPerComponent', 4);
```

This is intentionally separate from `main_spinodal_gnn_prep` / `downsample_graph_deterministic`.
It builds the graph from component-wise vertical midpoint pairing:
- each boundary connected component is handled independently
- boundary nodes are grouped in x-columns
- midpoint nodes are created from vertical pairs (`PairMode='extremes'` gives one vertical midpoint per selected column)
- a few representative midpoint nodes are kept per component
- representative nodes are connected only inside each component (MST), avoiding long global cross-connections

Distance attributes are stored per midpoint node (`VerticalGap`, `Distance`).

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
