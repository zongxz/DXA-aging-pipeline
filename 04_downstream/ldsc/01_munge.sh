# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# LDSC 第 1 步：GWAS 汇总统计 munge（对齐 HapMap3 位点）
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
LDSC=/path/to/shared_project/gwas/LDreg/ldsc
HM3=/path/to/shared_project/gwas/LDreg/data/eur_w_ld_chr/w_hm3.snplist
OUT=/path/to/project/GWAS/LDSC

$LDSC/munge_sumstats.py \
    --sumstats /path/to/project/GWAS/gwas_result/double_right.fastGWA \
    --out ${OUT}/double_right \
    --merge-alleles ${HM3} \
    --chunksize 500000

$LDSC/munge_sumstats.py \
    --sumstats /path/to/project/GWAS/gwas_result/single_right_knee.fastGWA \
    --out ${OUT}/single_right_knee \
    --merge-alleles ${HM3} \
    --chunksize 500000

$LDSC/munge_sumstats.py \
    --sumstats /path/to/project/GWAS/gwas_result/single_right_hip.fastGWA \
    --out ${OUT}/single_right_hip \
    --merge-alleles ${HM3} \
    --chunksize 500000