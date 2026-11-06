# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# plink2 逐染色体计算 PRS（GWAS 与非 GWAS 两套）
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
BASE_DIR="/path/to/project/GWAS/PRS_CS"
cd $BASE_DIR

for pheno in single_right_knee; do
    echo "===== Calculating PRS for ${pheno} ====="
    
    for chr in {1..22}; do
        plink2 \
        --bfile ${BASE_DIR}/geno_fixed/chr${chr} \
        --keep ${pheno}/GWAS.keep \
        --score ${pheno}/fixed/${pheno}_PRSCS_all_fixed.txt 2 4 6 cols=+scoresums \
        --out ${pheno}/score_gwas/chr${chr}_score
    done
done



for pheno in single_right_knee; do
    echo "===== Calculating PRS for ${pheno} ====="
    
    for chr in {1..22}; do
        plink2 \
        --bfile ${BASE_DIR}/geno_fixed/chr${chr} \
        --keep ${pheno}/nonGWAS.keep \
        --score ${pheno}/fixed/${pheno}_PRSCS_all_fixed.txt 2 4 6 cols=+scoresums \
        --out ${pheno}/score/chr${chr}_score
    done
done