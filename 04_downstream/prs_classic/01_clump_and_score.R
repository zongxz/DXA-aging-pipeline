# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# 基于 ieugwasr 的 LD clumping 与经典 PRS 计算
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
library(data.table)
library(dplyr)
library(ieugwasr)
library(genetics.binaRies)

setwd("/path/to/project/GWAS")

# 1. 读取 GWAS
knee <- fread("gwas_result/single_right_knee.txt")
hip <- fread("gwas_result/single_right_hip.txt")
double <- fread("gwas_result/double_right.txt")

# 2. 整理成 ld_clump 所需格式
knee <- knee %>% rename(rsid = SNP, pval = P)
hip <- hip %>% rename(rsid = SNP, pval = P)
double <- double %>% rename(rsid = SNP, pval = P)

# 3. 参数

bfile <- "/path/to/ref/1kg.v3/EUR"

phen_list <- c("single_right_knee", "single_right_hip", "double_right")

gwas_list <- list(
  single_right_knee = knee,
  single_right_hip  = hip,
  double_right      = double
)

cutoff_list <- c(0.05, 5e-8)

# 4. LD clumping

for (phen in phen_list) {
  gwas <- gwas_list[[phen]]
  
  cat("\n====================================\n")
  cat("Phenotype:", phen, "\n")
  cat("====================================\n")
  
  for (cutoff in cutoff_list) {
    cat("\nP <", cutoff, "\n")
    
    # P value 筛选
    sig <- gwas %>% filter(pval < cutoff)
    
    cat("SNPs before clumping:", nrow(sig), "\n")
    
    if (nrow(sig) == 0) {
      cat("No SNPs, skip.\n")
      next
    }
    
    # LD clumping
    clump_iv <- ld_clump(
      dat = sig,
      clump_kb = 500,
      clump_r2 = 0.1,
      clump_p = 1,
      plink_bin = genetics.binaRies::get_plink_binary(),
      bfile = bfile,
      pop = "EUR"
    )
    
    cat("SNPs after clumping:", nrow(clump_iv), "\n")
    
    # 检查重复
    cat("Duplicated SNPs:", sum(duplicated(clump_iv$rsid)), "\n")
    
    # SNP ID
    
    write.table(
      clump_iv$rsid,
      file = paste0("PRS/", phen, "_clump_", cutoff, ".txt"),
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
    
    # 完整 clumped GWAS
    write.table(
      clump_iv,
      file = paste0("PRS/", phen, "_gwas_clump_", cutoff, ".txt"),
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
  }
}