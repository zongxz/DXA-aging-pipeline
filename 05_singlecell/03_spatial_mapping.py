# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 05_singlecell
# cell2location 空间反卷积（Cell2location 模型）
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
import os
import gc
import torch
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp

from cell2location.models import Cell2location


# Paths

REF_FILE = ("/path/to/project/SingleCell/cell2location/reference_model/reference_posterior.h5ad")
SPATIAL_FILE = ("/path/to/project/SingleCell/cell2location/spatial_genes.h5ad")
OUT_DIR = ("/path/to/project/SingleCell/cell2location/spatial_model")

os.makedirs(OUT_DIR, exist_ok=True)


# GPU
print("=" * 70)
print("GPU")
print("=" * 70)

print("PyTorch:", torch.__version__)
print("CUDA:", torch.version.cuda)
print("CUDA available:", torch.cuda.is_available())

if not torch.cuda.is_available():
  raise RuntimeError("CUDA is not available.")

print("GPU:", torch.cuda.get_device_name(0))

torch.cuda.set_device(0)

# A800 Tensor Core optimization
torch.set_float32_matmul_precision("high")


# Load reference posterior
print("\n" + "=" * 70)
print("Loading reference posterior")
print("=" * 70)

adata_ref = sc.read_h5ad(REF_FILE)

print(adata_ref)

if "means_per_cluster_mu_fg" not in adata_ref.varm.keys():
  raise RuntimeError("means_per_cluster_mu_fg not found in reference posterior.")


# Extract reference signatures
inf_aver = adata_ref.varm["means_per_cluster_mu_fg"].copy()

print("\nReference signature matrix:")
print("shape:", inf_aver.shape)
print("type:", type(inf_aver))

print("\nCell-type columns:")
print(inf_aver.columns.tolist())


# Validate reference signatures
if inf_aver.isna().any().any():
  raise ValueError("Reference signatures contain NA.")

if (inf_aver.values < 0).any():
  raise ValueError("Reference signatures contain negative values.")


# Load spatial
print("\n" + "=" * 70)
print("Loading spatial data")
print("=" * 70)

adata_sp = sc.read_h5ad(SPATIAL_FILE)

print(adata_sp)

print("\nSpatial coordinates:",adata_sp.obsm["spatial"].shape)

if "spatial" not in adata_sp.obsm:
  raise RuntimeError("obsm['spatial'] not found.")

# Check spatial counts
X = adata_sp.X

if sp.issparse(X):
  values = X.data
else:
  values = np.asarray(X).ravel()

print("\nSpatial count range:")
print("min:", values.min())
print("max:", values.max())

if not np.all(values >= 0):
  raise ValueError("Spatial matrix contains negative values.")

if not np.all(values == np.round(values)):
  raise ValueError("Spatial matrix is not integer-valued raw counts.")

print("Raw-count check: PASS")


# Check genes

print("\n" + "=" * 70)
print("Checking genes")
print("=" * 70)

print("Reference genes:", adata_ref.n_vars)
print("Spatial genes:", adata_sp.n_vars)
print("Reference signature genes:", inf_aver.shape[0])

if not np.array_equal(adata_ref.var_names, inf_aver.index):
  raise RuntimeError("Reference var_names and inf_aver index differ.")

if not np.array_equal(adata_ref.var_names, adata_sp.var_names):
  raise RuntimeError("Reference and spatial genes are not identical or not in the same order.")

print("Gene matching: PASS")


# Cell2location setup

print("\n" + "=" * 70)
print("Cell2location setup")
print("=" * 70)

Cell2location.setup_anndata(adata=adata_sp)

# Create model
print("\nCreating Cell2location model...")

mod = Cell2location(
  adata_sp,
  cell_state_df=inf_aver,
  N_cells_per_location=30,
  detection_alpha=20
)

print("Model created.")


# Train
print("\n" + "=" * 70)
print("Training spatial model")
print("=" * 70)

mod.train(
  max_epochs=30000,
  batch_size=None,
  train_size=1.0
)

print("\nTraining finished.")


# Export posterior
print("\n" + "=" * 70)
print("Exporting posterior")
print("=" * 70)

adata_sp = mod.export_posterior(
  adata_sp,
  sample_kwargs={
    "num_samples": 1000,
    "batch_size": adata_sp.n_obs
  }
)

print("Posterior exported.")


# Save posterior
posterior_file = os.path.join(OUT_DIR, "spatial_posterior.h5ad")
adata_sp.write_h5ad(posterior_file)
print("\nPosterior saved:", posterior_file)


# Save model
model_dir = os.path.join(OUT_DIR,"model")
mod.save(model_dir, overwrite=True)
print("Model saved:", model_dir)

# Show posterior outputs
print("\n" + "=" * 70)
print("Posterior outputs")
print("=" * 70)

print("\nobsm:")
print(list(adata_sp.obsm.keys()))

print("\nuns:")
print(list(adata_sp.uns.keys()))

print("\nvarm:")
print(list(adata_sp.varm.keys()))

# Cleanup
del mod
del adata_ref
del adata_sp

gc.collect()

torch.cuda.empty_cache()

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)