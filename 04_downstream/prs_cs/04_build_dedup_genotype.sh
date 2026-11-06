# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# 构建去重基因型：重复 rsID 重命名为 rsID_A1_A2
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
GENO_PREFIX="/path/to/ref/info_pass" 
BASE_DIR="/path/to/project/GWAS/PRS_CS"


# 输出目录改到 BASE_DIR 下，避免污染原始数据 
OUT_DIR="${BASE_DIR}/geno_fixed" 
mkdir -p ${OUT_DIR}

for chr in {1..22}; do
    echo "===== Processing chr${chr} ====="
    
    # 1. 找出该染色体中重复的 rsID
    awk '{print $2}' ${GENO_PREFIX}/${chr}.bim | sort | uniq -d > geno_fixed/dup_chr${chr}.txt
    
    # 2. 将重复 rsID 重命名为 rsID_A1_A2 格式，生成新 bim
    awk '
    NR==FNR {dup[$1]=1; next}
    {
        if($2 in dup) $2=$2"_"$5"_"$6
        print
    }' OFS="\t" geno_fixed/dup_chr${chr}.txt ${GENO_PREFIX}/${chr}.bim > geno_fixed/chr${chr}.bim
    
    # 3. 复制配套的 bed 和 fam 文件
    cp ${GENO_PREFIX}/${chr}.bed geno_fixed/chr${chr}.bed
    cp ${GENO_PREFIX}/${chr}.fam geno_fixed/chr${chr}.fam
    
    # 4. 检查该染色体有多少重复SNP
    echo "chr${chr} duplicate SNPs: $(wc -l < geno_fixed/dup_chr${chr}.txt)"
done