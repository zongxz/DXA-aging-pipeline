#!/bin/bash
# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 03_gwas
# 基因型 QC 第 2 步：缺失率 / MAF / HWE / INFO 四重过滤
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------

out_folder=/path/to/project/GWAS
info_folder=/path/to/ref/ukb/ukb_imp_mfi/info_pass

for chr in {1..22}
do
    echo "=============================="
    echo "Processing chromosome ${chr}"
    echo "=============================="

    # 1. SNP missingness <= 5%
    plink2 \
      --bfile ${out_folder}/QC/${chr} \
      --geno 0.05 \
      --make-bed \
      --out ${out_folder}/info_pass/${chr}

    # 2. MAF >= 0.01
    # 重新计算一次，确保 missingness QC 后仍满足 MAF
    plink2 \
      --bfile ${out_folder}/info_pass/${chr} \
      --maf 0.01 \
      --make-bed \
      --out ${out_folder}/info_pass/${chr}

    # 3. HWE P >= 1e-6
    plink2 \
      --bfile ${out_folder}/info_pass/${chr} \
      --hwe 1e-6 \
      --make-bed \
      --out ${out_folder}/info_pass/${chr}

    # 4. INFO >= 0.8
    plink2 \
      --bfile ${out_folder}/info_pass/${chr} \
      --extract ${info_folder}/${chr}.txt \
      --make-bed \
      --out ${out_folder}/info_pass/${chr}

done

echo "All QC finished."