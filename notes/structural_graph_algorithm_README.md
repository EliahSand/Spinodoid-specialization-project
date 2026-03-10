# Structural Graph Extraction: Method Notes

This note explains why the current structural-graph pipeline works better, what was changed, and the exact algorithm from:

`dense reference graph -> skeleton graph -> structural graph`

Implementation source:

- `Matlab/gnn_prep_spinodal/src/graph/extract_structural_graph.m`

## 1) Why the current method works better

The current pipeline is better because it separates three jobs clearly:

1. **Topology extraction** (from dense to 1-pixel skeleton)  
2. **Model reduction** (compress degree-2 chains)  
3. **Local graph cleanup** (remove junction clutter and tiny artifacts)

Main improvements over earlier behavior:

1. **Topology-first reduction**  
   We do not directly downsample dense nodes. We first extract a topology-preserving skeleton, then reduce.

2. **Component-aware connectivity**  
   Skeleton edges are only added within the same connected component, and diagonal links are accepted only when supported by orthogonal neighbors. This prevents false bridges across nearby but disconnected spinodal regions.

3. **Chain-based compression instead of naive point decimation**  
   Nodes are not removed uniformly. Degree-2 chains are traced and compressed into structural edges while preserving junctions/endpoints.

4. **Post-pass contraction near junctions**  
   Very short edges and nearby junction clusters are merged in a controlled way. This removes tiny triangles and local node clutter.

5. **Geometry and provenance are preserved**  
   Each reduced edge keeps:
   - polyline coordinates along the original skeleton path,
   - mapping to original dense node ids,
   - physical path length (not just straight-line endpoint distance).

6. **Thickness is carried through**  
   Radius/thickness is computed from boundary distance, attached to nodes, and summarized per structural edge.

---

## 2) Data model (what each graph is)

Let the dense graph be

$$

G_d = (V_d, E_d), \quad |V_d| = N_d.

$$

Let the skeleton graph be

$$

G_s = (V_s, E_s), \quad |V_s| = N_s.

$$

Let the structural graph be

$$

G_r = (V_r, E_r), \quad |V_r| = N_r,\; N_r \ll N_s \le N_d.

$$

Node coordinates are in physical XY coordinates.

---

## 3) Dense graph -> raster occupancy

Dense nodes are binned to the native XY grid (clustered unique X and Y coordinates), producing:

- `occMask(r,c)`: occupied grid cells
- `nodeAt(r,c)`: dense node index at each occupied cell
- `componentAt(r,c)`: dense connected-component id
- `boundaryMask(r,c)`: from dense boundary nodes (fallback = geometric perimeter)

This gives a clean binary domain for morphology operations.

---

## 4) Skeleton extraction (Zhang-Suen thinning)

For each connected occupied component, we run **Zhang-Suen thinning** on the binary mask.

Credit:

- T. Y. Zhang and C. Y. Suen, *A Fast Parallel Algorithm for Thinning Digital Patterns*, Communications of the ACM, 1984.

### 4.1 Local neighborhood notation

For a foreground pixel $p_1$, neighbors $p_2,\dots,p_9$ follow the standard clockwise order.

Define:

$$

N(p_1) = \sum_{i=2}^{9} p_i

$$

$$

S(p_1) = \text{number of } 0 \to 1 \text{ transitions in } (p_2,p_3,\dots,p_9,p_2).

$$

### 4.2 Deletion rules (two sub-iterations)

A pixel is deleted only if all hold:

$$

2 \le N(p_1) \le 6,\qquad S(p_1)=1

$$

and:

- **Sub-iteration 1**:
  $$

  \neg(p_2 p_4 p_6),\qquad \neg(p_4 p_6 p_8)
  
$$
- **Sub-iteration 2**:
  $$

  \neg(p_2 p_4 p_8),\qquad \neg(p_2 p_6 p_8)
  
$$

Repeat until no deletions occur.

Result: a one-node-thick centerline skeleton that preserves topology.

### 4.3 Small-island safeguard

If thinning removes a tiny component too aggressively, we enforce at least `MinIslandNodes` by keeping the most interior pixel (max distance from boundary) in that component.

---

## 5) Thickness/radius computation

On occupied cells:

1. Compute distance to boundary seeds in physical XY units.
2. Define:

$$

r(x,y) = d_{\partial \Omega}(x,y), \qquad t(x,y) = 2r(x,y).

$$

These are attached to:

- dense nodes (`node_radius`, `node_thickness`),
- skeleton nodes (inherited from mapped dense nodes),
- structural nodes (inherited from mapped skeleton nodes).

---

## 6) Skeleton graph construction

Skeleton pixels become nodes in $G_s$.  
Edges are added between neighboring skeleton pixels using local directions.

Important guard against false links:

- A diagonal link is allowed only if at least one orthogonal support cell exists in the same component.
- This blocks corner-touch shortcuts across disconnected spinodal regions.

Node role from degree:

- degree $=0$: island
- degree $=1$: endpoint
- degree $=2$: chain/waypoint
- degree $>2$: junction

---

## 7) Skeleton -> structural compression

Compression traces maximal paths through degree-2 chains and replaces each chain segment by one reduced edge.

Anchors are always retained at endpoints/junctions. Extra waypoints are controlled by `DetailLevel`.

### 7.1 Open chain detail rule

For a traced path with $n$ nodes:

$$

n_{\text{internal}} = n-2

$$

$$

n_{\text{extra}} = \text{round}(\text{DetailLevel} \cdot n_{\text{internal}})

$$

clipped to $[0, n_{\text{internal}}]$.

Anchors are:

- first node,
- `n_extra` evenly spaced internal positions,
- last node.

So:

- `DetailLevel = 0`: only chain endpoints kept (coarsest)
- `DetailLevel = 1`: all internal chain nodes kept (finest, near-skeleton resolution)

### 7.2 Closed loop detail rule

For a cycle with $n$ nodes:

$$

n_{\text{keep}} = \min\left(n,\; \max\left(3,\; 1+\text{round}(\text{DetailLevel}(n-1))\right)\right)

$$

At least 3 anchors are kept to preserve a loop representation.

### 7.3 What each reduced edge stores

For each structural edge $e$:

- endpoint indices in structural graph,
- `edge_length`: polyline arc length
  $$

  L_e = \sum_k \|p_{k+1}-p_k\|_2
  
$$
- `edge_polyline`: ordered skeleton coordinates,
- `edge_full_node_paths`: ordered dense-node ids along the path.

This is why reduction does not lose geometric traceability.

---

## 8) Cleanup pass (why clutter is reduced)

After initial compression, we run a minimal local contraction pass.

### 8.1 Characteristic spacing

Estimate grid spacing $h$ (median axis spacing; fallback to median edge length).

Set:

$$

\ell_{\min} = 1.5h,\qquad r_{\text{merge}} = 1.5h

$$

with a small fixed iteration cap.

### 8.2 Operations per iteration

1. **Short-edge contraction**: merge edge endpoints if edge length $< \ell_{\min}$ and eligibility checks pass.
2. **Junction-cluster merge**: merge nearby high-degree nodes (same component, distance $< r_{\text{merge}}$).

### 8.3 Representative node choice when merging

Choose representative in this priority:

1. highest degree,
2. highest thickness (if available),
3. medoid in XY (minimum total distance).

Then:

- rewire incident edges to representative,
- remove self-loops,
- deduplicate parallel edges (keep strongest/longest),
- preserve mapping sets of merged original indices.

This specifically removes tiny triangles, clustered junction duplicates, and very short redundant edges.

---

## 9) Structural edge thickness summaries

Using each edge path in dense indices, compute:

$$

t_{\min},\; t_{\text{mean}},\; t_{\max},\; t_{\text{bottleneck}}=t_{\min},\; \Delta t=t_{\max}-t_{\min}.

$$

These are stored per structural edge.

---

## 10) DetailLevel: practical interpretation

`DetailLevel` controls only how many intermediate waypoints are retained on degree-2 chains/cycles.

- Low values (for example 0.1 to 0.25):
  - smaller graph,
  - faster training/inference,
  - less local geometric detail.

- Mid values (about 0.3 to 0.5):
  - good compromise for most runs.

- High values (0.7 to 1.0):
  - close to skeleton resolution,
  - larger graph, more local detail.

Important: topology comes from thinning + chain tracing; `DetailLevel` mainly controls resolution, not connectivity class.

### 10.1 DetailLevel mathematics (exact implementation)

Let

- $d \in [0,1]$ be `DetailLevel`.

#### Open chain

For a traced open chain with $n$ skeleton nodes:

$$
m = n - 2
$$

where $m$ is the number of internal degree-2 nodes.

The number of extra internal anchors kept is:

$$
a = \min\!\big(m,\;\mathrm{round}(d\,m)\big).
$$

Anchors kept on that chain:

$$
N_{\mathrm{anchors}} = 2 + a
$$

(first endpoint + last endpoint + $a$ internal anchors).

Reduced structural edges created from that chain:

$$
E_{\mathrm{chain}} = N_{\mathrm{anchors}} - 1 = a + 1.
$$

Boundary cases:

- $d=0 \Rightarrow a=0$: only chain endpoints retained (maximum compression),
- $d=1 \Rightarrow a=m$: all internal chain nodes retained (no chain compression).

#### Closed cycle

For a cycle with $n$ nodes:

$$
k = \min\!\left(n,\; \max\!\left(3,\; 1+\mathrm{round}(d\,(n-1))\right)\right),
$$

where $k$ is the number of cycle anchors retained.

So the algorithm always keeps at least 3 loop anchors and at most all $n$ cycle nodes.

#### Global intuition

If $K$ is the number of key nodes (junctions/endpoints/islands) and $C$ is the number of degree-2 chain nodes, a useful approximation is:

$$
N_{\mathrm{struct}} \approx K + d\,C
$$

(with small deviations due to loop minimum-anchor constraints and cleanup merges).

### 10.2 Why this works (simple summary)

In simple terms, the method works better because:

1. It simplifies **after** extracting a topology-preserving centerline (skeleton), not directly on dense mesh nodes.
2. It always preserves key shape logic (junctions/endpoints), and only compresses redundant degree-2 path segments.
3. `DetailLevel` gives a predictable coarse-to-fine control over chain resolution.
4. A cleanup pass removes local clutter (very short edges, clustered junction duplicates) without changing large-scale connectivity.
5. Geometry/physics traceability is preserved through per-edge polyline paths, path lengths, node-id mappings, and thickness summaries.

---

## 11) Why this is stable for ML graph inputs

The resulting structural graph is more ML-friendly because:

1. Topology is preserved from skeletonization.
2. Redundant local clutter is removed by deterministic cleanup.
3. Edge lengths and paths are physically meaningful.
4. Thickness information is explicitly available at node and edge level.
5. Mapping back to skeleton/dense nodes is retained for auditability and feature backfilling.
