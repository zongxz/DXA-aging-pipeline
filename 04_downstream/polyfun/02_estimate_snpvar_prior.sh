# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# 第 2 步：--compute-h2-L2，估计全基因组 SNP 先验方差
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
cd /path/to/project/GWAS/PolyFun

python software/polyfun/polyfun.py \
    --compute-h2-L2 \
    --allow-missing \
    --output-prefix ldscores/double_right \
    --sumstats munge/double_right.parquet \
    --ref-ld-chr data/baselineLF2.2.UKB/baselineLF2.2.UKB. \
    --w-ld-chr data/baselineLF2.2.UKB/weights.UKB. \
    2>&1 | tee logs/double_right_compute_h2_L2.log