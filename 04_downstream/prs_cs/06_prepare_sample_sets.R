# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# 划分 GWAS 人群 / 非 GWAS 人群（.keep 文件）
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
# ============================================
# R script: prepare_ids.R
# ============================================

library(data.table)

# ========== 基础路径配置 ==========
fam_path <- "/path/to/ref/info_pass/1.fam"  # 任选一个染色体的 fam 即可
out_dir  <- "/path/to/project/GWAS/PRS_CS"

# ========== 三个表型的文件路径配置 ==========
GWAS_PATH <- '/path/to/project/data/gwas'
pheno_config <- list(
  double_right = list(
    phen   = file.path(GWAS_PATH, 'double_right.txt'),
    covar  = file.path(GWAS_PATH, 'covar.txt'),
    qcovar = file.path(GWAS_PATH, 'qcovar.txt')    
  ),
  single_right_knee = list(
    phen   = file.path(GWAS_PATH, 'single_right_knee.txt'), 
    covar  = file.path(GWAS_PATH, 'covar.txt'),
    qcovar = file.path(GWAS_PATH, 'qcovar.txt')      
  )
  # single_right_hip = list(
  #   phen   = file.path(GWAS_PATH, 'single_right_hip.txt'),
  #   covar  = file.path(GWAS_PATH, 'covar.txt'),
  #   qcovar = file.path(GWAS_PATH, 'qcovar.txt')    
  # )
)

# 读取 fam（所有染色体共用）
fam <- fread(fam_path, header = FALSE, select = 1:2, col.names = c("FID", "IID"))
fam_ids <- paste(fam$FID, fam$IID)
cat("Total genotyped samples:", length(fam_ids), "\n\n")

# 循环三个表型
for(pheno_name in names(pheno_config)){
  
  cfg <- pheno_config[[pheno_name]]
  cat("===== Processing:", pheno_name, "=====\n")
  
  # 读取表型、协变量
  phen   <- fread(cfg$phen,   header = FALSE, na.strings = c("NA", "-9"), col.names = c("FID", "IID", "phenotype"))
  covar  <- fread(cfg$covar,  header = FALSE, na.strings = c("NA", "-9"))
  qcovar <- fread(cfg$qcovar, header = FALSE, na.strings = c("NA", "-9"))
  
  # 提取无缺失样本
  phen_ids   <- phen[complete.cases(phen),   .(FID, IID)]
  covar_ids  <- covar[complete.cases(covar), .(V1, V2)]
  qcovar_ids <- qcovar[complete.cases(qcovar), .(V1, V2)]
  
  setnames(covar_ids,  c("V1", "V2"), c("FID", "IID"))
  setnames(qcovar_ids, c("V1", "V2"), c("FID", "IID"))
  
  # 取交集：表型 + 协变量 + 数量协变量
  analysis_ids <- Reduce(intersect, list(
    paste(phen_ids$FID, phen_ids$IID),
    paste(covar_ids$FID, covar_ids$IID),
    paste(qcovar_ids$FID, qcovar_ids$IID)
  ))
  
  # 再与基因型样本取交集 → 最终 GWAS 样本
  gwas_ids <- intersect(fam_ids, analysis_ids)
  cat("  GWAS samples:", length(gwas_ids), "\n")
  
  # 保存 GWAS 人群 ID
  gwas_df <- data.frame(
    FID = sub(" .*", "", gwas_ids),
    IID = sub(".* ", "", gwas_ids)
  )
  fwrite(gwas_df, file.path(out_dir, pheno_name, "GWAS.txt"), 
         sep = "\t", col.names = TRUE)
  fwrite(gwas_df, file.path(out_dir, pheno_name, "GWAS.keep"), 
       sep = "\t", col.names = FALSE)
  # 非 GWAS 人群 = 基因型中有，但不在该表型 GWAS 中
  non_gwas_ids <- setdiff(fam_ids, gwas_ids)
  cat("  non-GWAS samples:", length(non_gwas_ids), "\n")
  
  non_gwas_df <- data.frame(
    FID = sub(" .*", "", non_gwas_ids),
    IID = sub(".* ", "", non_gwas_ids)
  )
  
  # 保存到各自目录下，避免覆盖
  fwrite(non_gwas_df, file.path(out_dir, pheno_name, "nonGWAS.keep"), 
         sep = "\t", col.names = FALSE)
  
  all <- data.frame(
    FID = sub(" .*", "", fam_ids),
    IID = sub(".* ", "", fam_ids)
  )
  
  cat("  ALL samples:", nrow(all), "\n")
  fwrite(all, file.path(out_dir, pheno_name, "ALL.keep"), 
         sep = "\t", col.names = FALSE)
  
  cat("  Saved to:", file.path(out_dir, pheno_name), "\n\n")
}

  cat("All phenotypes completed!\n")