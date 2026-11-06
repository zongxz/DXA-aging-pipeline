# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# 合并 22 条染色体的后验效应文件
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
cd /path/to/project/GWAS/PRS_CS/PRSCS_output

for trait in single_right_knee
do
    echo "===== ${trait} ====="

    for chr in {1..22}
    do
        cat ${trait}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt
    done > ${trait}_PRSCS_all_chr.txt

done