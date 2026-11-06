# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 03_gwas
# GCTA fastGWA-mlm 关联分析（含稀疏 GRM 与协变量）
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
input_folder=/path/to/project/GWAS
out_folder=/path/to/project/GWAS


covar=/path/to/project/data/gwas/covar.txt
qcovar=/path/to/project/data/gwas/qcovar.txt
sp_grm=${out_folder}/sp_grm/grm

phen=/path/to/project/data/gwas/single_right_knee.txt


geno_chrs_auto=/path/to/project/data/gwas/geno_chrs_auto.txt
head $geno_chrs_auto



/path/to/tools/gcta/gcta64 --mbfile ${geno_chrs_auto} \
--grm-sparse ${sp_grm} \
--fastGWA-mlm \
--covar ${covar} \
--qcovar ${qcovar} \
--pheno ${phen} \
--nofilter \
--thread-num 32 \
--out ${out_folder}/gwas_result/fastGWA