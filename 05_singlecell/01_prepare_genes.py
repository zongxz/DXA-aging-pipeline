# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 05_singlecell
# 读入单细胞参考与空间转录组，做基因过滤
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
import os
import scanpy as sc
import numpy as np

REF_FILE = "/path/to/project/SingleCell/data/GSE169396_scRNA.h5ad"
SP_FILE = "/path/to/project/SingleCell/data/GSE284089_spatial.h5ad"

OUT_DIR = "/path/to/project/SingleCell/cell2location"
os.makedirs(OUT_DIR, exist_ok=True)

REF_OUT = os.path.join(OUT_DIR, "reference_genes.h5ad")

SP_OUT = os.path.join(OUT_DIR, "spatial_genes.h5ad")

# 1. Load
print("Loading reference...")
adata_ref = sc.read_h5ad(REF_FILE)

print("Reference:", adata_ref)

print("\nLoading spatial...")
adata_sp = sc.read_h5ad(SP_FILE)

print("Spatial:", adata_sp)

# 2. Check genes
print("\n========== Gene information ==========")
print("Reference genes:", adata_ref.n_vars)
print("Spatial genes:", adata_sp.n_vars)
print("Reference duplicated genes:", adata_ref.var_names.duplicated().sum())
print("Spatial duplicated genes:", adata_sp.var_names.duplicated().sum())

# 3. Make gene names unique if necessary

if adata_ref.var_names.duplicated().any():
    adata_ref.var_names_make_unique()

if adata_sp.var_names.duplicated().any():
    adata_sp.var_names_make_unique()

# 4. Find common genes

common_genes = adata_ref.var_names[adata_ref.var_names.isin(adata_sp.var_names)]

print("\n========== Gene intersection ==========")
print("Common genes:", len(common_genes))
print("Reference:", adata_ref.n_vars)

print("Spatial:", adata_sp.n_vars)
print("Reference coverage:", len(common_genes) / adata_ref.n_vars)
print("Spatial coverage:", len(common_genes) / adata_sp.n_vars)

# 5. Keep reference gene order
common_genes = adata_ref.var_names[adata_ref.var_names.isin(adata_sp.var_names)]

adata_ref = adata_ref[:, common_genes].copy()

adata_sp = adata_sp[:, common_genes].copy()

# 6. Force identical gene order

adata_sp = adata_sp[:, adata_ref.var_names].copy()

assert np.array_equal(adata_ref.var_names,adata_sp.var_names)

# 7. Save
adata_ref.write_h5ad(REF_OUT)
adata_sp.write_h5ad(SP_OUT)

# 8. Final check
print("\n========== Final ==========")
print("Reference:",adata_ref)

print("Spatial:",adata_sp)

print("Same genes:",np.array_equal(adata_ref.var_names,adata_sp.var_names))

print("\nSaved:")
print(REF_OUT)
print(SP_OUT)
