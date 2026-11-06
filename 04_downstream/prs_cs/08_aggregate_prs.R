# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# 合并 22 条染色体的 PRS 得分
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
library(data.table)

base_dir <- "/path/to/project/GWAS/PRS_CS"
phenotypes <- c("single_right_knee", "double_right")

for(pheno in phenotypes){
    cat("===== Merging GWAS", pheno, "=====\n")
    
    prs <- NULL
    
    for(chr in 1:22){
        file_path <- file.path(base_dir, pheno, "score_gwas", paste0("chr", chr, "_score.sscore"))
        
        if(!file.exists(file_path)){
            warning(paste("Missing file:", file_path))
            next
        }
        
        x <- fread(file_path, header=TRUE)
        x <- x[, .(IID, SCORE1_SUM)]
        setnames(x, "SCORE1_SUM", paste0("chr", chr))
        
        if(is.null(prs)){
            prs <- x
        } else {
            prs <- merge(prs, x, by="IID", all=TRUE)
        }
    }
    
    # 计算总PRS（22条染色体之和）
    chr_cols <- paste0("chr", 1:22)
    prs[, PRS := rowSums(.SD, na.rm=TRUE), .SDcols=chr_cols]
    
    # 保存完整结果（含各染色体分数）
    fwrite(prs, file.path(base_dir, pheno, paste0(pheno, "_PRS_GWAS_full.csv")))
    
    # 保存简洁结果（仅IID和总PRS）
    fwrite(prs[, .(IID, PRS)], file.path(base_dir, pheno, paste0(pheno, "_PRS_GWAS.csv")))
    
    cat(pheno, "done. N =", nrow(prs), "\n")

    cat("===== Merging non-GWAS", pheno, "=====\n")
    prs <- NULL
    
    for(chr in 1:22){
        file_path <- file.path(base_dir, pheno, "score", paste0("chr", chr, "_score.sscore"))
        
        if(!file.exists(file_path)){
            warning(paste("Missing file:", file_path))
            next
        }
        
        x <- fread(file_path, header=TRUE)
        x <- x[, .(IID, SCORE1_SUM)]
        setnames(x, "SCORE1_SUM", paste0("chr", chr))
        
        if(is.null(prs)){
            prs <- x
        } else {
            prs <- merge(prs, x, by="IID", all=TRUE)
        }
    }
    
    # 计算总PRS（22条染色体之和）
    chr_cols <- paste0("chr", 1:22)
    prs[, PRS := rowSums(.SD, na.rm=TRUE), .SDcols=chr_cols]
    
    # 保存完整结果（含各染色体分数）
    fwrite(prs, file.path(base_dir, pheno, paste0(pheno, "_PRS_non-GWAS_full.csv")))
    
    # 保存简洁结果（仅IID和总PRS）
    fwrite(prs[, .(IID, PRS)], file.path(base_dir, pheno, paste0(pheno, "_PRS_non-GWAS.csv")))
    
    cat(pheno, "done. N =", nrow(prs), "\n")
}

cat("All phenotypes completed!\n")