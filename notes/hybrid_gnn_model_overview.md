# Hybrid GNN Model Overview

This note explains how the hybrid model in `Matlab/GNN/models/Main_hybrid_vision.m`
works. The model combines two different deep learning modules:

- a graph neural network branch for structural connectivity
- a convolutional neural network branch for dense raster geometry

The final output is a vector of PCA coefficients that reconstructs the midpoint
out-of-plane displacement profile `u3(s)`.

## High-Level Idea

The hybrid model uses two complementary views of the same geometry.

The graph branch sees the structure as nodes and edges. This is useful for
learning from connectivity, local neighborhoods, graph topology, node radius,
and boundary information.

The dense branch sees the structure as a 2D raster image. This is useful for
learning spatial patterns, global layout, holes, boundaries, and image-like
geometric features.

Both branches produce fixed-length feature vectors. These vectors are
concatenated together, along with global design parameters, and passed through
fully connected fusion layers. The fused representation is decoded into the PCA
target.

In compact form:

```text
graph geometry  -> GNN branch -> graph feature vector
dense geometry  -> CNN branch -> dense feature vector
global params   -> normalized features

[graph features; dense features; global features]
        -> fusion MLP
        -> decoder
        -> predicted PCA coefficients
        -> reconstructed u3(s)
```

## Model Modes

The script supports three modes:

```matlab
modelMode = 'hybrid';
modelMode = 'dense_only';
modelMode = 'graph_only';
```

In `hybrid` mode, both the graph branch and dense raster branch are active.

In `graph_only` mode, the dense branch is replaced by zeros, so the model uses
only graph features plus global parameters.

In `dense_only` mode, the graph branch is replaced by zeros, so the model uses
only raster features plus global parameters.

This is useful for ablation studies, because the same fusion and decoder code
can compare whether the graph branch, dense branch, or combination performs
best.

## Inputs

Each sample has three input groups.

### 1. Structural Graph Input

The graph input is loaded as:

```text
X_cell{g}  : F x Ng node feature matrix
ei_cell{g} : 2 x Eg edge index matrix
N_vec(g)   : number of valid nodes
```

For the hybrid model, the node features are:

```text
[x, y, radius, boundary]
```

where:

- `x, y` are node coordinates
- `radius` describes the local structural radius/thickness feature
- `boundary` is a binary flag marking boundary nodes

The graph edges define which nodes are connected.

Before training, graphs are padded to the same maximum number of nodes:

```text
X_pad   : 4 x maxN x nSamples
nodeMask: 1 x maxN x nSamples
```

The mask is needed because different graphs can have different numbers of
nodes. Padded dummy nodes must not contribute to graph pooling.

### 2. Dense Raster Input

The dense input is a rasterized image-like representation of the geometry:

```text
Dense: H x W x 2 x nSamples
```

The two stored channels are:

```text
channel 1: occupancy
channel 2: boundary
```

Before entering the CNN, two coordinate channels are added:

```text
channel 3: x_grid
channel 4: y_grid
```

So the CNN input becomes:

```text
Dense_input: H x W x 4 x nSamples
```

The coordinate channels help the CNN know where a feature appears in the
domain. Without coordinate channels, a convolutional network is mostly
translation-equivariant, which is useful for images but not always ideal for
geometry where absolute position matters.

### 3. Global Design Parameters

The model also receives two global parameters:

```text
thickness ratio
angle
```

Thickness ratio is standardized using training-set statistics:

```matlab
tr_std = (tr_ratio - tr_mean) / tr_std
```

Angle is not used directly as degrees. It is encoded as:

```matlab
sin(2 * angle)
cos(2 * angle)
```

where `angle` is first converted to radians.

The final global feature vector is therefore:

```text
[standardized thickness ratio,
 sin(2 angle),
 cos(2 angle)]
```

The factor `2` represents the symmetry/periodicity of the sampled angle. This
also avoids problems where angles that are physically close could look far apart
as raw numbers.

## Graph Branch

The graph branch is implemented in `hybrid_forward.m` as `graph_branch`.

Input shape:

```text
X: 4 x maxN x B
```

where:

- `4` is the number of node features
- `maxN` is the padded node count
- `B` is the batch size

The current model configuration uses:

```text
K = 4
hiddenDim = 96
```

### Node Embedding

Each node feature vector is first projected into a hidden representation:

```matlab
u = relu(W_embedding * X + b_embedding)
```

Shape:

```text
u: hiddenDim x maxN x B
```

With the current configuration:

```text
u: 96 x maxN x B
```

The node mask is applied after embedding so padded nodes are zeroed:

```matlab
u = u .* mask
```

### Normalized Adjacency

Each graph has a normalized adjacency matrix:

```matlab
A_hat = D^(-1/2) (A + I) D^(-1/2)
```

where:

- `A` is the graph adjacency matrix
- `I` adds self-loops
- `D` is the degree matrix of `A + I`

The self-loop means each node keeps its own information during message passing.
The symmetric normalization keeps feature magnitudes more stable across nodes
with different degrees.

### Message Passing Layers

For each graph layer `i = 1,...,K`, the model computes:

```matlab
v = W_i * u + b_i
Vagg = A_hat * v
u = relu(u + Vagg)
u = u .* mask
```

The important part is:

```text
new node state = old node state + aggregated neighbor message
```

This is a residual graph update. The residual connection helps preserve
previous information while adding neighborhood information.

After one layer, a node has information from its immediate neighbors. After
four layers, a node can contain information from a wider four-hop neighborhood.

### Graph Pooling

The model needs one fixed-length vector per graph, not one vector per node.
Therefore it pools node embeddings into graph-level features.

After the initial embedding and after each graph layer, the model computes:

```matlab
hMean = mean of valid node embeddings
hMax  = max of valid node embeddings
hPool = [hMean; hMax]
```

The model stores pooled features from:

```text
embedding output
after graph layer 1
after graph layer 2
after graph layer 3
after graph layer 4
```

Since there are `K + 1` pooled states, and each state contributes both mean and
max pooling, the graph feature dimension is:

```text
graphDim = 2 * hiddenDim * (K + 1)
```

With the current values:

```text
graphDim = 2 * 96 * (4 + 1) = 960
```

So the graph branch outputs:

```text
hGraph: 960 x 1 x B
```

## Dense CNN Branch

The dense branch is implemented in `hybrid_forward.m` as `dense_branch`.

Input shape:

```text
Dense: H x W x 4 x B
```

For the current raster setup:

```text
Dense: 128 x 128 x 4 x B
```

The channels are:

```text
occupancy, boundary, x_grid, y_grid
```

The current CNN configuration is:

```text
cnnChannels = [32, 64, 128]
```

The CNN applies three convolutional layers:

```text
Conv 1: 5 x 5, stride 2,  32 channels, ReLU
Conv 2: 3 x 3, stride 2,  64 channels, ReLU
Conv 3: 3 x 3, stride 2, 128 channels, ReLU
```

In MATLAB code:

```matlab
z = dlconv(Dense, CNN1.W, CNN1.b, 'Padding', 'same', 'Stride', 2);
z = relu(z);
z = dlconv(z, CNN2.W, CNN2.b, 'Padding', 'same', 'Stride', 2);
z = relu(z);
z = dlconv(z, CNN3.W, CNN3.b, 'Padding', 'same', 'Stride', 2);
z = relu(z);
```

With a 128 x 128 input and stride 2 in each layer, the approximate spatial
sizes are:

```text
128 x 128 -> 64 x 64 -> 32 x 32 -> 16 x 16
```

The final CNN feature map has 128 channels.

### CNN Pooling

The CNN output is converted to a fixed-length vector using both global mean
pooling and global max pooling:

```matlab
hMean = mean over spatial dimensions
hMax  = max over spatial dimensions
hDense = [hMean; hMax]
```

The dense branch dimension is:

```text
cnnDim = 2 * finalCNNChannels
```

With the current model:

```text
cnnDim = 2 * 128 = 256
```

So the dense branch outputs:

```text
hDense: 256 x 1 x B
```

## How the Two Deep Learning Modules Are Combined

The GNN and CNN are not combined by mixing their internal layers. They are
combined at the feature-vector level.

Each branch first processes its own input type independently:

```text
graph input -> GNN -> hGraph
raster input -> CNN -> hDense
```

The result from each branch is a fixed-length learned representation of the same
sample.

Then the model concatenates:

```matlab
h = cat(1, hGraph, hDense, G_global);
```

This means the full fused vector is:

```text
h = [graph representation;
     dense raster representation;
     global parameter representation]
```

With the current dimensions:

```text
hGraph  = 960 features
hDense  = 256 features
G_global = 3 features
```

Therefore:

```text
fusion input dimension = 960 + 256 + 3 = 1219
```

This concatenation is the key operation that combines the two deep learning
modules. After concatenation, the following fully connected layers learn how to
weight, mix, and use information from both branches.

Conceptually:

- the GNN says: "Here is what I learned from structural connectivity"
- the CNN says: "Here is what I learned from the rasterized geometry"
- the global vector says: "Here are the design parameters"
- the fusion MLP learns how to combine all three sources into one prediction

This is a common multimodal learning strategy. Each input modality gets a
specialized encoder, and the encoded features are fused before prediction.

## Fusion MLP

After concatenation, the fused vector is passed through two dense layers:

```matlab
h = relu(Fusion1.W * h + Fusion1.b);
h = relu(Fusion2.W * h + Fusion2.b);
```

The current fusion dimension is:

```text
fusionDim = 192
```

So the fusion path is:

```text
1219 -> 192 -> 192
```

Dropout is applied during training:

```text
dropoutRate = 0.30
```

Dropout is applied to:

- graph branch output
- dense branch output
- first fusion hidden vector

Dropout is disabled during validation and prediction.

## Decoder

The final decoder is a linear layer:

```matlab
Zhat = Decoder.W * h + Decoder.b
```

The output dimension is:

```text
nComp = 8
```

So the model outputs:

```text
Zhat: 8 x 1 x B
```

These are predicted standardized PCA coefficients.

## Target and Reconstruction

The neural network does not directly output the full 128-point `u3(s)` curve.
It outputs PCA coefficients.

The PCA target has:

```text
8 coefficients per sample
```

The training target is standardized using training-set statistics:

```matlab
Z_target = (Z - Z_mean) ./ Z_std
```

After prediction, coefficients are unstandardized:

```matlab
Z_pred = Zhat_std .* Z_std + Z_mean
```

Then the physical displacement curve is reconstructed:

```matlab
U3_pred = Z_pred * coeff' + u3_mean
```

where:

- `coeff` is the PCA basis
- `u3_mean` is the training-set mean profile
- `U3_pred` is the reconstructed 128-point displacement profile

## Loss Function

The model is trained in standardized PCA space.

The loss is weighted mean squared error:

```matlab
err = (Zhat - Z_target).^2;
loss = mean(sum(w .* err, 1), 'all');
```

The weights `w` are based on explained variance from PCA:

```matlab
loss_w = explained(1:nComp);
loss_w = loss_w / mean(loss_w);
```

This gives more importance to PCA components that explain more variance in the
deformation profiles.

## Training Objective and Early Stopping

The model is optimized with Adam.

Current main training settings:

```text
batchSize    = 16
maxEpochs    = 180
lr0          = 7e-4
decay_rate   = 0.985
weight_decay = 1e-4
valFreq      = 5
patience     = 8
```

Although the loss is computed in PCA space, the best checkpoint is selected
using validation `u3` R2:

```text
best model = model with highest validation u3 R2
```

This matters because the final engineering quantity of interest is the physical
displacement curve, not just the PCA coefficients.

## Evaluation Metrics

The model is evaluated in PCA space and physical displacement space.

PCA-space metrics:

```text
per-PC R2
weighted PCA R2
unweighted PCA R2
```

Physical-space metrics:

```text
u3 MSE
u3 R2
```

The most interpretable metric is `u3 R2`, because it measures the reconstructed
displacement profile directly.

## Why This Hybrid Design Makes Sense

The graph branch and dense branch capture different information.

The GNN is good for:

- structural connectivity
- local graph neighborhoods
- node-level radius and boundary information
- topology-sensitive features

The CNN is good for:

- spatial patterns in the full geometry
- image-like layout information
- broad global shape cues
- occupancy and boundary distributions

The fusion MLP learns how to combine these learned summaries. It does not assume
in advance whether graph information or dense information is more important.
That weighting is learned from the data during training.

In other words, the model uses:

```text
specialized encoders first, shared prediction head second
```

This is the central mechanism that allows two different deep learning modules to
work together in one model.

## Exact Current Dimension Summary

Using the current `Main_hybrid_vision.m` configuration:

```text
Graph input node features: 4
Graph hidden dimension:   96
Graph layers:             4
Graph pooled states:      K + 1 = 5
Graph pooling:            mean + max
Graph output dimension:   2 * 96 * 5 = 960

Dense input channels:     4
CNN channels:             [32, 64, 128]
Dense pooling:            mean + max
Dense output dimension:   2 * 128 = 256

Global features:          3
Fusion input dimension:   960 + 256 + 3 = 1219
Fusion hidden dimension:  192
Decoder output dimension: 8 PCA coefficients
```

Full architecture:

```text
Graph branch:
4 x maxN
-> node embedding to 96
-> 4 residual graph message-passing layers
-> mean/max pooling at 5 depths
-> 960 graph features

Dense branch:
128 x 128 x 4
-> Conv 5x5, 32, stride 2
-> Conv 3x3, 64, stride 2
-> Conv 3x3, 128, stride 2
-> global mean/max pooling
-> 256 dense features

Global branch:
[tr_ratio, angle]
-> [standardized tr_ratio, sin(2 angle), cos(2 angle)]
-> 3 global features

Fusion:
[960 graph; 256 dense; 3 global]
-> 1219 fused features
-> dense 192 + ReLU
-> dense 192 + ReLU
-> linear decoder
-> 8 standardized PCA coefficients
```

