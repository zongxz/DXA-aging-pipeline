# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 05_singlecell
# 提取细胞丰度与后验结果
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
import os
import scanpy as sc
import pandas as pd
import numpy as np

# ============================================================
# Paths
# ============================================================

input_file = (
    "/path/to/project/"
    "SingleCell/cell2location/spatial_model/spatial_posterior.h5ad"
)

outdir = (
    "/path/to/project/"
    "SingleCell/cell2location/spatial_model/results"
)

os.makedirs(outdir, exist_ok=True)

# ============================================================
# Load posterior
# ============================================================

print("=" * 70)
print("Loading spatial posterior")
print("=" * 70)

adata = sc.read_h5ad(input_file)

print(adata)

# ============================================================
# 1. Cell abundance
# ============================================================

abundance = adata.obsm[
    "means_cell_abundance_w_sf"
].copy()

# Clean cell type names
abundance.columns = [
    x.split("means_per_cluster_mu_fg_")[-1]
    for x in abundance.columns
]

print("\nCell abundance:")
print("Shape:", abundance.shape)
print("Cell types:")
print(abundance.columns.tolist())

# ============================================================
# 2. Uncertainty
# ============================================================

q05 = adata.obsm[
    "q05_cell_abundance_w_sf"
].copy()

q95 = adata.obsm[
    "q95_cell_abundance_w_sf"
].copy()

std = adata.obsm[
    "stds_cell_abundance_w_sf"
].copy()

q05.columns = abundance.columns
q95.columns = abundance.columns
std.columns = abundance.columns

# ============================================================
# 3. Cell proportion
# ============================================================

row_sum = abundance.sum(axis=1)

proportion = abundance.div(
    row_sum.replace(0, np.nan),
    axis=0
)

print("\nCell proportion:")
print("Shape:", proportion.shape)

print("\nProportion row sums:")
print(proportion.sum(axis=1).describe())

# ============================================================
# 4. Put cleaned results into AnnData
# ============================================================

adata.obsm["cell_abundance"] = abundance
adata.obsm["cell_abundance_q05"] = q05
adata.obsm["cell_abundance_q95"] = q95
adata.obsm["cell_abundance_std"] = std
adata.obsm["cell_proportion"] = proportion

# ============================================================
# 5. Export CSV
# ============================================================

abundance_file = os.path.join(
    outdir,
    "cell_abundance.csv"
)

proportion_file = os.path.join(
    outdir,
    "cell_proportion.csv"
)

q05_file = os.path.join(
    outdir,
    "cell_abundance_q05.csv"
)

q95_file = os.path.join(
    outdir,
    "cell_abundance_q95.csv"
)

std_file = os.path.join(
    outdir,
    "cell_abundance_std.csv"
)

abundance.to_csv(abundance_file)
proportion.to_csv(proportion_file)
q05.to_csv(q05_file)
q95.to_csv(q95_file)
std.to_csv(std_file)

# ============================================================
# 6. Save updated h5ad
# ============================================================

output_h5ad = os.path.join(
    outdir,
    "spatial_posterior_with_proportion.h5ad"
)

adata.write(output_h5ad)

# ============================================================
# 7. Final validation
# ============================================================

print("\n" + "=" * 70)
print("FINAL VALIDATION")
print("=" * 70)

print("Spots:", adata.n_obs)
print("Genes:", adata.n_vars)
print("Cell types:", abundance.shape[1])

print("\nCell types:")
print(abundance.columns.tolist())

print("\nAbundance:")
print(abundance.head())

print("\nProportion:")
print(proportion.head())

print("\nNaN abundance:", abundance.isna().sum().sum())
print("Negative abundance:", (abundance < 0).sum().sum())

print(
    "Proportion row-sum deviation:",
    np.abs(proportion.sum(axis=1) - 1).max()
)

print("\nSaved files:")
print(abundance_file)
print(proportion_file)
print(q05_file)
print(q95_file)
print(std_file)
print(output_h5ad)

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)