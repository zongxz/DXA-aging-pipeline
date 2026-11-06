# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 05_singlecell
# cell2location 参考模型（RegressionModel）训练，得到细胞类型表达签名
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
import os
import gc
import torch
import scanpy as sc
from cell2location.models import RegressionModel

# Configuration

REF_FILE = ("/path/to/project/SingleCell/cell2location/reference_genes.h5ad")

OUT_DIR = ("/path/to/project/SingleCell/cell2location/reference_model")

os.makedirs(OUT_DIR, exist_ok=True)

# GPU
print("PyTorch:", torch.__version__)
print("CUDA:", torch.version.cuda)
print("CUDA available:", torch.cuda.is_available())

if torch.cuda.is_available():
    print("GPU:", torch.cuda.get_device_name(0))
    torch.cuda.set_device(0)

# Load reference

adata_ref = sc.read_h5ad(REF_FILE)

print("\nReference:")
print(adata_ref)
print("\nCell types:")
print(adata_ref.obs["cell_type"].value_counts())
print("\nSamples:")
print(adata_ref.obs["sample"].value_counts())

# Check raw counts
import scipy.sparse as sp
import numpy as np

X = adata_ref.X

if sp.issparse(X):
    values = X.data
else:
    values = np.asarray(X).ravel()
print("\nCount range:", values.min(), values.max())

assert np.all(values >= 0)
assert np.all(values == np.round(values))

print("Raw count check: PASS")

# Setup AnnData
RegressionModel.setup_anndata(
    adata=adata_ref,
    labels_key="cell_type",
    batch_key="sample"
)

# Create model
mod = RegressionModel(adata_ref)

print("\nModel created.")

# Train
mod.train(
    max_epochs=250,
    batch_size=2500,
    train_size=1.0
)

# Export posterior

print("\nExporting posterior...")

adata_ref = mod.export_posterior(
    adata_ref,
    sample_kwargs={
        "num_samples": 1000,
        "batch_size": 2500
    }
)

# Save posterior
posterior_file = os.path.join(OUT_DIR, "reference_posterior.h5ad")

adata_ref.write_h5ad(posterior_file)

print("\nPosterior saved:",posterior_file)

# Save model

model_dir = os.path.join(OUT_DIR, "model")

mod.save(model_dir, overwrite=True)

print("Model saved:", model_dir)

# Cleanup
del mod
del adata_ref

gc.collect()

if torch.cuda.is_available():
    torch.cuda.empty_cache()

print("\nReference training finished.")
