# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# 创建各表型输出目录
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
BASE_DIR="/path/to/project/GWAS/PRS_CS"
cd $BASE_DIR

# 创建各表型输出目录
# score_gwas / score 分别存放 GWAS 人群与非 GWAS 人群的打分结果，缺一 plink2 会直接报错
for pheno in double_right single_right_knee single_right_hip; do
    mkdir -p ${pheno}/score ${pheno}/score_gwas ${pheno}/fixed
done

# 基因型去重后的存放目录（三个表型共用）
mkdir -p geno_fixed