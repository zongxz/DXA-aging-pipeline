# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# 逐染色体运行 PRS-CS，推断后验效应量
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
for trait in single_right_hip double_right single_right_knee
do

    echo "=============================================="
    echo "Trait ${trait} started at $(date)"
    echo "=============================================="

    for chr in {1..22}
    do
        echo "===== ${trait}: Chromosome ${chr} started at $(date) ====="

        python /path/to/project/GWAS/PRS_CS/PRScs/PRScs.py \
            --ref_dir=/path/to/ref/PRSCS_ref/ldblk_1kg_eur \
            --bim_prefix=/path/to/project/GWAS/info_pass/${chr} \
            --sst_file=/path/to/project/GWAS/gwas_result/${trait}.prscs.txt \
            --n_gwas=48850 \
            --chrom=${chr} \
            --out_dir=/path/to/project/GWAS/PRS_CS/PRSCS_output/${trait}

        if [ $? -eq 0 ]; then
            echo "===== ${trait}: Chromosome ${chr} finished ====="
        else
            echo "===== ${trait}: Chromosome ${chr} FAILED ====="
        fi
    done

    echo "=============================================="
    echo "Trait ${trait} finished at $(date)"
    echo "=============================================="

done