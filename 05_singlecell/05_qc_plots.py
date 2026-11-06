# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 05_singlecell
# 反卷积结果 QC 图
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# ============================================================
# Paths
# ============================================================

input_file = (
    "/path/to/project/"
    "SingleCell/cell2location/spatial_model/spatial_posterior.h5ad"
)

outdir = (
    "/path/to/project/"
    "SingleCell/cell2location/spatial_model/spatial_QC"
)

os.makedirs(outdir, exist_ok=True)

# ============================================================
# Load
# ============================================================

adata = sc.read_h5ad(input_file)

print(adata)

# ============================================================
# Extract cell2location abundance
# ============================================================

abundance = adata.obsm["means_cell_abundance_w_sf"].copy()

# Clean names
abundance.columns = [
    x.split("means_per_cluster_mu_fg_")[-1]
    for x in abundance.columns
]

print("\nCell types:")
print(abundance.columns.tolist())

# ============================================================
# Spatial coordinates
# ============================================================

coords = pd.DataFrame(
    adata.obsm["spatial"],
    index=adata.obs_names,
    columns=["x", "y"]
)

print("\nCoordinates:")
print(coords.head())

# ============================================================
# Combine
# ============================================================

plot_data = coords.join(abundance)

# ============================================================
# Save abundance
# ============================================================

plot_data.to_csv(
    os.path.join(outdir, "cell2location_abundance_with_coordinates.csv")
)

# ============================================================
# Plot each cell type
# ============================================================

for celltype in abundance.columns:

    fig, ax = plt.subplots(figsize=(8, 8))

    sc = ax.scatter(
        plot_data["x"],
        plot_data["y"],
        c=plot_data[celltype],
        s=12,
        cmap="viridis",
        linewidths=0
    )

    ax.set_title(
        f"cell2location abundance: {celltype}",
        fontsize=14
    )

    ax.set_xlabel("x")
    ax.set_ylabel("y")

    # Spatial image coordinates usually have y increasing downward
    ax.invert_yaxis()

    ax.set_aspect("equal")

    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("Estimated abundance")

    plt.tight_layout()

    outfile = os.path.join(
        outdir,
        f"abundance_{celltype}.png"
    )

    plt.savefig(
        outfile,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print("Saved:", outfile)

# ============================================================
# Summary statistics
# ============================================================

summary = abundance.describe().T

summary["zero_fraction"] = (
    (abundance == 0).sum(axis=0) /
    abundance.shape[0]
).values

summary.to_csv(
    os.path.join(outdir, "abundance_summary.csv")
)

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print(summary)

print("\nDONE")