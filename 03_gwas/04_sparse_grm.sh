# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 03_gwas
# GCTA 稀疏 GRM 分块并行构建与合并
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
out_folder=/path/to/project/GWAS
hapmap_snp=${out_folder}/hapmap/hapmap.prune.in
head $hapmap_snp

geno_chrs_auto=/path/to/project/data/gwas/geno_chrs_auto.txt
geno_chrs=${geno_chrs_auto}

for i in {1..180}
do
  echo "i=$i"
  /path/to/tools/gcta/gcta64 --mbfile ${geno_chrs} --extract ${hapmap_snp} \
  --make-grm-part 180 $i --sparse-cutoff 0.05 --thread-num 64 \
  --out ${out_folder}/sp_grm/sp_grm
done  

# Merge all the parts together (Linux, Mac)
cat ${out_folder}/sp_grm/sp_grm.part_180_*grm.id > ${out_folder}/sp_grm/grm.grm.id
cat ${out_folder}/sp_grm/sp_grm.part_180_*grm.sp > ${out_folder}/sp_grm/grm.grm.sp

wc -l ${out_folder}/sp_grm/grm.grm.sp
wc -l ${out_folder}/sp_grm/grm.grm.id