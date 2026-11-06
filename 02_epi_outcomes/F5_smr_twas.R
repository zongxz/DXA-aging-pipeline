# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 02_epi_outcomes
# SMR 与 TWAS 整合分析
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
# SMR ====
library(data.table)

folder_path <- 'result/F5_SMR/'
files <- list.files(path = folder_path, pattern = "\\.tsv$")

prefixes <- unique(sub("\\d+\\.tsv$", "", files))

for (p in prefixes) {
  # 匹配该前缀的所有文件
  group_files <- grep(paste0("^", p, "\\d+\\.tsv$"), files, value = TRUE)
  group_files <- paste0(folder_path, group_files)
  merged <- rbindlist(lapply(group_files, fread))
  
  out_name <- paste0(folder_path, p, ".tsv")
  fwrite(merged, out_name, sep = "\t")
  cat("已合并生成:", out_name, "\n")
}


# TWAS ====
library(data.table)

folder_path <- 'result/F5_TWAS/result/'
files <- list.files(path = folder_path, pattern = "\\.csv$")

prefixes <- c("double_right", "single_right_hip", "single_right_knee")

# 循环合并每个部位
for (pre in prefixes) {
  target_files <- paste0(folder_path, grep(paste0("^", pre, "_"), files, value = TRUE))
  
  combined <- map_dfr(target_files, function(f) {
    df <- read.csv(f)
    df$tissue <- gsub(paste0("^", pre, "_"), "", f)  # 提取组织名
    df$tissue <- gsub("\\.csv$", "", df$tissue)      # 去掉后缀
    return(df)
  })
  
  # 保存合并后的文件
  write.csv(combined, paste0('result/F5_TWAS/', pre, "_TWAS.csv"), row.names = FALSE)
  
  cat("已合并:", pre, "->", length(target_files), "个组织,", nrow(combined), "行\n")
}




# SMR FDR ====
folder_path <- 'result/F5_SMR/'
files <- list.files(path = folder_path, pattern = "\\.tsv$")
# filter(p_HEIDI > 0.05) 

for (file in files) {
  df <- fread(file.path(folder_path, file))
  names(df)[1] <- 'Gene'

  df <- df %>%
    filter(qtl_name == 'eQTL_GTEx_Whole_Blood') %>%
    filter(p_HEIDI > 0.05) %>%
    mutate(
      FDR_SMR = p.adjust(p_SMR, method = "BH"),
      zscore = b_SMR / se_SMR,
      CI_lower = b_SMR - qnorm(0.975) * se_SMR,
      CI_upper = b_SMR + qnorm(0.975) * se_SMR,
      tissue = gsub('eQTL_GTEx_', '', qtl_name)
    ) %>%
    filter(FDR_SMR < 0.05) %>%
    dplyr::select(tissue, Gene, probeID, gene_id, 
                  b_SMR, se_SMR, zscore,CI_lower,CI_upper,
                  p_SMR, p_SMR_multi, p_HEIDI, FDR_SMR)
  
  write.csv(df, gsub('.tsv', '_FDR_Whole_Blood.csv', file.path(folder_path, file)), row.names = F)
    
}

# TWAS FDR ====
folder_path <- 'result/F5_TWAS'
files <- list.files(path = folder_path, pattern = "\\.csv$")
# filter(p_HEIDI > 0.05) 

for (file in files) {
  df <- fread(file.path(folder_path, file))

  df <- df %>%
    filter(grepl('Whole_Blood', tissue)) %>%
    mutate(
      effect_size = as.numeric(effect_size),
      se = abs(effect_size / zscore),
      CI_lower = effect_size - qnorm(0.975) * se,
      CI_upper = effect_size + qnorm(0.975) * se,
      FDR_TWAS = p.adjust(pvalue, method = 'BH'),
      tissue = gsub('result/F5_TWAS/result/', '', tissue)
    ) %>%
    mutate(tissue = gsub(gsub('TWAS.csv', '', file), '', tissue)) %>%
    filter(FDR_TWAS < 0.05) %>%
    dplyr::select(tissue, gene, gene_name, zscore, se,
                  effect_size, pvalue, FDR_TWAS,CI_lower,CI_upper)
  
  print(table(df$tissue))
  write.csv(df, gsub('.csv', '_FDR_Whole_Blood.csv', file.path(folder_path, file)), row.names = F)
}
