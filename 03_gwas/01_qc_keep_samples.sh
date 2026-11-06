# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 03_gwas
# 基因型 QC 第 1 步：按分析样本清单保留个体
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
input_folder=/path/to/raw_genotype
out_folder=/path/to/project/GWAS
keep_ind=/path/to/project/data/gwas/gwas_id.txt

mkdir -p ${out_folder}/QC

for chr in {1..22}
do
    echo "Processing chromosome ${chr}"

    plink2 \
      --bfile ${input_folder}/${chr}_v0 \
      --keep ${keep_ind} \
      --maf 0.01 \
      --make-bed \
      --out ${out_folder}/QC/${chr}
done