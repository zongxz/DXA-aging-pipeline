# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# PolyFun 第 1 步：汇总统计 munge 为 parquet
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
cd /path/to/project/GWAS/PolyFun

mkdir -p munge prior ldscores logs

#python software/polyfun/munge_polyfun_sumstats.py \
#    --sumstats /path/to/project/GWAS/gwas_result/single_right_knee.txt.gz \
#    --out munge/single_right_knee.parquet \
#    --min-info 0 \
#    --min-maf 0.001

python software/polyfun/munge_polyfun_sumstats.py \
    --sumstats /path/to/project/GWAS/gwas_result/double_right.txt.gz \
    --out munge/double_right.parquet \
    --min-info 0 \
    --min-maf 0.001