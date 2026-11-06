# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# 同步重命名权重文件中的重复 rsID
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
BASE_DIR="/path/to/project/GWAS/PRS_CS"
cd $BASE_DIR/PRSCS_output

# 先合并所有染色体的重复rsID映射关系（只需生成一次）
cat $BASE_DIR/geno_fixed/dup_chr*.txt | sort | uniq > $BASE_DIR/all_dup_rsid.txt

for pheno in single_right_knee; do
    echo "===== Fixing ${pheno} ====="
    
    # 生成 rename 映射表: 原始rsID -> rsID_A1_A2
    # 从 PRSCS_all_chr.txt 中提取重复SNP的 A1,A2 信息
    awk '
    NR==FNR {dup[$1]=1; next}
    dup[$2] {print $2, $2"_"$4"_"$5}
    ' $BASE_DIR/all_dup_rsid.txt ${pheno}_PRSCS_all_chr.txt > $BASE_DIR/${pheno}/fixed/rename_map.txt
    
    # 对 PRSCS_all_chr.txt 进行重命名
    awk '
    NR==FNR {map[$1]=$2; next}
    {
        if($2 in map) $2=map[$2]
        print
    }' OFS="\t" $BASE_DIR/${pheno}/fixed/rename_map.txt ${pheno}_PRSCS_all_chr.txt \
    > $BASE_DIR/${pheno}/fixed/${pheno}_PRSCS_all_fixed.txt
    
    echo "${pheno} fixed file: $BASE_DIR/${pheno}/fixed/${pheno}_PRSCS_all_fixed.txt"
done