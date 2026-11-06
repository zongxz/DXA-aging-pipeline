# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 03_gwas
# HapMap3 位点提取与 LD 剪枝（用于稀疏 GRM）
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
out_folder=/path/to/project/GWAS
hapmap=/path/to/ref/ukb/hapmap/hapmap3_r1_b36_fwd.CEU.qc.poly.recode

plink \
--file ${hapmap} \
--maf 0.01 \
--indep-pairwise 1000 100 0.9 \
--out ${out_folder}/hapmap

