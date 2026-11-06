# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# TWAS 第 2 步：逐组织运行 S-PrediXcan
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
for trait in double_right single_right_hip single_right_knee
do

    echo "============================================================"
    echo "Trait: ${trait}"
    echo "Started at: $(date)"
    echo "============================================================"

    cat /path/to/project/GWAS/TWAS/tissue.txt | while read tissue;
    do

        echo "------------------------------------------------------------"
        echo "Trait  : ${trait}"
        echo "Tissue : ${tissue}"
        echo "Start  : $(date)"
        echo "------------------------------------------------------------"

        python /path/to/other_project/twas/software/SPrediXcan.py \
            --gwas_file /path/to/project/GWAS/TWAS/${trait}.txt.gz \
            --snp_column snp \
            --effect_allele_column effect_allele \
            --non_effect_allele_column non_effect_allele \
            --zscore_column zscore \
            --model_db_path /path/to/other_project/twas/data/models/eqtl/mashr/mashr_${tissue}.db \
            --covariance /path/to/other_project/twas/data/models/eqtl/mashr/mashr_${tissue}.txt.gz \
            --keep_non_rsid \
            --additional_output \
            --model_db_snp_key varID \
            --throw \
            --output_file /path/to/project/GWAS/TWAS/result/${trait}_${tissue}.csv

        if [ $? -eq 0 ]; then
            echo "Finished successfully: ${trait} - ${tissue}"
        else
            echo "FAILED: ${trait} - ${tissue}"
        fi

    done

    echo "============================================================"
    echo "Trait ${trait} finished at $(date)"
    echo "============================================================"

done