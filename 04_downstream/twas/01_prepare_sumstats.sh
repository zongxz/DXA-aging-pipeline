# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# TWAS 第 1 步：GWAS 位点映射到预测模型 SNP 并生成 betas
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
for trait in double_right single_right_hip single_right_knee
do

    echo "===== ${trait} started at $(date) ====="

    gwas=/path/to/project/GWAS/gwas_result/${trait}.fastGWA

    python /path/to/other_project/twas/software/M03_betas.py \
        --snp_map_file /path/to/other_project/twas/data/map_snp150_hg38.txt.gz \
        --gwas_file ${gwas} \
        --snp_column SNP \
        --non_effect_allele_column A2 \
        --effect_allele_column A1 \
        --beta_column BETA \
        --pvalue_column P \
        --keep_non_rsid \
        --throw \
        --output /path/to/project/GWAS/TWAS/${trait}.txt.gz

    if [ $? -eq 0 ]; then
        echo "===== ${trait} finished successfully at $(date) ====="
    else
        echo "===== ${trait} FAILED at $(date) ====="
    fi

done