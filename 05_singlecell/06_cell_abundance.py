#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 05_singlecell
# 细胞丰度汇总与可视化
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# 1. 输入文件
# ============================================================

abundance_file = (
    "/path/to/project/"
    "SingleCell/cell2location/spatial_model/results/"
    "cell_abundance.csv"
)

coord_file = (
    "/path/to/project/"
    "SingleCell/cell2location/spatial_model/spatial_QC/"
    "cell2location_abundance_with_coordinates.csv"
)

# 输出目录
outdir = (
    "/path/to/project/"
    "SingleCell/cell2location/spatial_model/results/"
    "dominant_celltype"
)

os.makedirs(outdir, exist_ok=True)


# ============================================================
# 2. 读取数据
# ============================================================

print("Loading abundance...")
abundance = pd.read_csv(abundance_file, index_col=0)

print("Abundance shape:", abundance.shape)
print("Cell types:")
print(abundance.columns.tolist())

print("\nLoading coordinates...")
coord = pd.read_csv(coord_file, index_col=0)

print("Coordinate shape:", coord.shape)
print("Coordinate columns:")
print(coord.columns.tolist())


# ============================================================
# 3. 统一 spot ID
# ============================================================

common = abundance.index.intersection(coord.index)

print("\nNumber of common spots:", len(common))

abundance = abundance.loc[common].copy()
coord = coord.loc[common].copy()


# ============================================================
# 4. 确保 cell2location abundance 全部为数值
# ============================================================

abundance = abundance.apply(pd.to_numeric, errors="coerce")

if abundance.isna().any().any():
    print("WARNING: NA values detected.")
    print(abundance.isna().sum())
    raise ValueError("Abundance contains NA values.")

if (abundance < 0).any().any():
    raise ValueError("Abundance contains negative values.")


# ============================================================
# 5. 每个 spot 找最高丰度 cell type
# ============================================================

# 第一高
dominant_celltype = abundance.idxmax(axis=1)
dominant_abundance = abundance.max(axis=1)

# 排序后得到第一高和第二高
sorted_abundance = np.sort(abundance.values, axis=1)

dominant_abundance_array = sorted_abundance[:, -1]
second_abundance_array = sorted_abundance[:, -2]

# 第二高
second_celltype = abundance.apply(
    lambda x: x.nlargest(2).index[1],
    axis=1
)

# 第一高 / 所有细胞类型总丰度
abundance_sum = abundance.sum(axis=1)

dominant_proportion = (
    dominant_abundance / abundance_sum
)

# 第一高与第二高之间的相对差异
dominance_margin = (
    (dominant_abundance - second_abundance_array)
    / dominant_abundance
)


# ============================================================
# 6. 构建结果表
# ============================================================

result = pd.DataFrame({
    "x": coord["x"],
    "y": coord["y"],
    "dominant_celltype": dominant_celltype,
    "dominant_abundance": dominant_abundance,
    "second_celltype": second_celltype,
    "second_abundance": second_abundance_array,
    "dominant_proportion": dominant_proportion,
    "dominance_margin": dominance_margin
})

# 保存
result_file = os.path.join(
    outdir,
    "dominant_celltype_per_spot.csv"
)

result.to_csv(result_file)

print("\nSaved:")
print(result_file)


# ============================================================
# 7. 每种 dominant cell type 数量
# ============================================================

celltype_counts = (
    result["dominant_celltype"]
    .value_counts()
    .rename("n_spots")
)

celltype_percent = (
    result["dominant_celltype"]
    .value_counts(normalize=True)
    .mul(100)
    .rename("percent")
)

summary = pd.concat(
    [celltype_counts, celltype_percent],
    axis=1
)

summary_file = os.path.join(
    outdir,
    "dominant_celltype_summary.csv"
)

summary.to_csv(summary_file)

print("\n============================================================")
print("Dominant cell type summary")
print("============================================================")
print(summary)

print("\nSaved:")
print(summary_file)


# ============================================================
# 8. dominant proportion 的统计
# ============================================================

prop_summary = result["dominant_proportion"].describe()

print("\n============================================================")
print("Dominant proportion")
print("============================================================")
print(prop_summary)

prop_file = os.path.join(
    outdir,
    "dominant_proportion_summary.csv"
)

prop_summary.to_csv(prop_file)


# ============================================================
# 9. dominance margin 的统计
# ============================================================

margin_summary = result["dominance_margin"].describe()

print("\n============================================================")
print("Dominance margin")
print("============================================================")
print(margin_summary)

margin_file = os.path.join(
    outdir,
    "dominance_margin_summary.csv"
)

margin_summary.to_csv(margin_file)


# ============================================================
# 10. 绘制 dominant cell type 空间分布
# ============================================================

plt.figure(figsize=(10, 9))

for celltype in abundance.columns:

    sub = result[
        result["dominant_celltype"] == celltype
    ]

    if len(sub) == 0:
        continue

    plt.scatter(
        sub["x"],
        sub["y"],
        s=25,
        alpha=0.85,
        label=f"{celltype} (n={len(sub)})"
    )

plt.gca().invert_yaxis()

plt.xlabel("x")
plt.ylabel("y")

plt.title(
    "Dominant cell type by cell2location abundance"
)

plt.legend(
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    frameon=False
)

plt.tight_layout()

plot_file = os.path.join(
    outdir,
    "dominant_celltype_spatial.png"
)

plt.savefig(
    plot_file,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("\nSaved:")
print(plot_file)


# ============================================================
# 11. 绘制 dominant proportion 空间图
# ============================================================

plt.figure(figsize=(9, 9))

sc = plt.scatter(
    result["x"],
    result["y"],
    c=result["dominant_proportion"],
    s=25,
    alpha=0.9
)

plt.gca().invert_yaxis()

plt.xlabel("x")
plt.ylabel("y")

plt.title(
    "Dominant cell type proportion"
)

plt.colorbar(
    sc,
    label="Dominant proportion"
)

plt.tight_layout()

plot_file = os.path.join(
    outdir,
    "dominant_proportion_spatial.png"
)

plt.savefig(
    plot_file,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Saved:")
print(plot_file)


# ============================================================
# 12. 绘制 dominance margin 空间图
# ============================================================

plt.figure(figsize=(9, 9))

sc = plt.scatter(
    result["x"],
    result["y"],
    c=result["dominance_margin"],
    s=25,
    alpha=0.9
)

plt.gca().invert_yaxis()

plt.xlabel("x")
plt.ylabel("y")

plt.title(
    "Dominance margin: top1 vs top2"
)

plt.colorbar(
    sc,
    label="(Top1 - Top2) / Top1"
)

plt.tight_layout()

plot_file = os.path.join(
    outdir,
    "dominance_margin_spatial.png"
)

plt.savefig(
    plot_file,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Saved:")
print(plot_file)


# ============================================================
# 13. 按 dominance proportion 分级
# ============================================================

def classify_confidence(x):

    if x >= 0.50:
        return "High"
    elif x >= 0.33:
        return "Moderate"
    else:
        return "Mixed"


result["dominance_confidence"] = (
    result["dominant_proportion"]
    .apply(classify_confidence)
)


confidence_summary = (
    result["dominance_confidence"]
    .value_counts()
    .rename("n_spots")
    .to_frame()
)

confidence_summary["percent"] = (
    confidence_summary["n_spots"]
    / len(result)
    * 100
)

confidence_file = os.path.join(
    outdir,
    "dominance_confidence_summary.csv"
)

confidence_summary.to_csv(confidence_file)

print("\n============================================================")
print("Dominance confidence")
print("============================================================")
print(confidence_summary)

print("\nSaved:")
print(confidence_file)


# ============================================================
# 14. 输出最终结果
# ============================================================

final_file = os.path.join(
    outdir,
    "dominant_celltype_per_spot_with_confidence.csv"
)

result.to_csv(final_file)

print("\n============================================================")
print("DONE")
print("============================================================")
print("Final result:")
print(final_file)

print("\nOutput directory:")
print(outdir)

