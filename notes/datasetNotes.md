# Dataset Notes

## Group-wise validation R2

The subgroup validation R2 values were computed after training from saved predictions, not during the MATLAB training loop. The latest hybrid run stores `valMetrics.U3_pred` in `Matlab/GNN/models/best/20260428_230347_hybrid/final_metrics.mat`. These predictions were matched to the original validation samples using `val_idx` from `Matlab/GNN/data/dataset/targets/pca_targets.mat`, then compared against `U3_mat` from `Matlab/GNN/data/dataset/targets/u3_targets.mat`.

For each subgroup, the same deformation-space R2 formula was recomputed only on that subset:

```text
R2 = 1 - sum((U3_pred - U3_true)^2) / sum((U3_true - mean(U3_true))^2)
```

Latest validation breakdown:

```text
By thickness ratio:
tr=0.50: 0.831
tr=1.00: 0.864
tr=2.00: 0.779

By angle band:
0-30 deg:  0.760
31-60 deg: 0.738
61-90 deg: 0.910
```

This suggests the current hybrid model generalizes worst for `tr=2.00` and for lower/mid angle ranges.