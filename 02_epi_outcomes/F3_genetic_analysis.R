# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 02_epi_outcomes
# 遗传分析：药物靶点、GWAS Catalog、FUMA、经典 PRS、LDpred2
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------



# 骨龄加速曼哈顿图 ====

library(data.table)
library(ggplot2)
library(patchwork)

# 文件
files <- c(double = "/path/to/epi/data/GWAS_scPagwas/double_right.txt", 
           hip    = "/path/to/epi/data/GWAS_scPagwas/single_right_hip.txt", 
           knee   = "/path/to/epi/data/GWAS_scPagwas/single_right_knee.txt")

files <- c(hip_SHom  = "/path/to/epi/data/GWAS_scPagwas/data_hip.txt.gz", 
           knee_SHom = "/path/to/epi/data/GWAS_scPagwas/data_knee.txt.gz", 
           NGWAMA    = "/path/to/epi/data/GWAS_scPagwas/knee_hip.N_weighted_GWAMA.results.txt.gz")

NGWAMA <- fread("/path/to/epi/data/GWAS_scPagwas/knee_hip.N_weighted_GWAMA.results.txt.gz")
names(NGWAMA)

NGWAMA <- NGWAMA[,c('CHR', 'BP','PVAL')]
names(NGWAMA) <- c('CHR', 'POS','P')

SHom <- fread("/path/to/epi/data/GWAS_scPagwas/data_hip.txt.gz")
names(SHom)

SHom <- SHom[,c('CHR', 'POS','P')]

# 读取数据
gwas_data <- list()

for (nm in names(files)) {
  cat("Reading:", nm, "\n")
  x <- fread(files[[nm]])
  
  cat("Rows:", nrow(x), "\n")
  cat("Columns:", paste(names(x), collapse = ", "), "\n")
  
  x$CHR <- as.integer(x$CHR)
  x$POS <- as.numeric(x$POS)
  x$P   <- as.numeric(x$P)
  
  # 只保留 chr1-22
  x <- x[!is.na(x$CHR) & x$CHR >= 1 & x$CHR <= 22 & 
           !is.na(x$POS) & !is.na(x$P) & x$P > 0 & x$P <= 1, ]
  
  setorder(x, CHR, POS)
  gwas_data[[nm]] <- x
  cat("After filtering:", nrow(x), "\n\n")
}


nm = 'NGWAMA'
gwas_data[[nm]] <- NGWAMA
# 获取三个 GWAS 中每条染色体的最大位置
chr_len <- data.frame(CHR = 1:22, chr_len = 0)

for (chr in 1:22) {
  max_positions <- c()
  
  for (nm in names(gwas_data)) {
    x <- gwas_data[[nm]]
    tmp <- x$POS[x$CHR == chr]
    
    if (length(tmp) > 0) {
      max_positions <- c(max_positions, max(tmp, na.rm = TRUE))
    }
  }
  
  if (length(max_positions) > 0) {
    chr_len$chr_len[chr] <- max(max_positions)
  }
}


# 4. 计算 genome-wide offset
# 加入染色体间 gap
gap <- 10e6

chr_len$offset <- 0

for (i in 2:22) {
  chr_len$offset[i] <- chr_len$offset[i - 1] +  chr_len$chr_len[i - 1] + gap
}

# 染色体中心
chr_len$center <- chr_len$offset + chr_len$chr_len / 2

# BPcum
for (nm in names(gwas_data)) {
  x <- gwas_data[[nm]]
  x$BPcum <- x$POS + chr_len$offset[x$CHR]
  gwas_data[[nm]] <- x
}

# 检查
print(chr_len)
sapply(gwas_data, function(x) max(x$BPcum))

# 绘图代码
make_manhattan <- function(x, show_x = FALSE) {
  color_palette <- c(
    "#DC143C", "#1E90FF", "#05004e", "#007b43", "#f8b500",
    "#8B008B", "#dbd0e6", "#FF69B4", "#ff847c", "#a0a003",
    "#08d9d6", "#ba52ed", "#323232", "#8B4513", "#00A6A6",
    "#6A5ACD", "#E67E22", "#2E8B57", "#C71585", "#708090",
    "#B8860B", "#4169E1"
  )
  
  p <- ggplot(x, aes(x = BPcum, y = -log10(P))) +
    geom_point(aes(color = factor(CHR)), size = 0.35, alpha = 0.6) +
    scale_color_manual(values = color_palette, guide = "none") +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", linewidth = 0.4) +
    scale_x_continuous(breaks = chr_len$center, labels = chr_len$CHR, expand = c(0.005, 0)) +
    labs(title = "", x = NULL, y = expression(bold('-') ~ bold(log)[bold("10")] ~ bold('(') ~ bolditalic('P') ~ bold(')'))) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0, face = "bold", family = "Times New Roman"),
      axis.title = element_text(size = 16, family = "Times New Roman", face = "bold"),
      axis.text.y = element_text(size = 14, family = "Times New Roman", face = "bold"),
      axis.text.x = if (show_x) {
        element_text(size = 14, family = "Times New Roman", face = "bold")
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x) {
        element_line()
      } else {
        element_blank()
      },
      plot.margin = margin(3, 8, 3, 8)
    )
}

p_knee <- make_manhattan(gwas_data$knee)
p_hip <- make_manhattan(gwas_data$hip)
p_double <- make_manhattan(gwas_data$double, show_x = TRUE)

final_plot <-p_hip / p_knee / p_double 

ggsave('result/F3_FUMA/manhattan.png', final_plot, width = 10, height = 10, dpi=300)

p_knee_SHom <- make_manhattan(gwas_data$knee_SHom, show_x = TRUE)

ggsave('result/F3_FUMA/manhattan_knee_SHom.png', p_knee_SHom, width = 10, height = 3, dpi=300)

p_NGWAMA <- make_manhattan(gwas_data$NGWAMA, show_x = TRUE)

ggsave('result/F3_FUMA/manhattan_NGWAMA.png', p_NGWAMA, width = 10, height = 3, dpi=300)
# N =====
# 韦恩图 ====
# 三组共定位方式基因与MAGMA基因

traits <- c('double_right', 'single_right_hip', 'single_right_knee')

for (trait in traits) {
  folder_path <- '/path/to/epi/result/F3_FUMA/'
  
  pos <- fread(paste0(folder_path, trait, '/POS/genes.txt'))
  eqtl <- fread(paste0(folder_path, trait, '/eQTL/genes.txt'))
  chr_3d <- fread(paste0(folder_path, trait, '/CHR3D/genes.txt'))
  magma <- fread(paste0(folder_path, trait, '/POS/magma.genes.out'))
  magma <- magma %>%
    filter(P <= 0.05 / nrow(magma))
  
  library(ggvenn)
  unique_gene <- unique(c(pos$ensg, eqtl$ensg, chr_3d$ensg, magma$GENE))
  common_genes <- Reduce(intersect, list(pos$ensg, eqtl$ensg, chr_3d$ensg, magma$GENE))
  venn_df <- data.frame(gene_ensg = unique_gene,
                        Positional = unique_gene %in% pos$ensg,
                        eQTL = unique_gene %in% eqtl$ensg,
                        Chromatin = unique_gene %in% chr_3d$ensg,
                        MAGMA = unique_gene %in% magma$GENE)
  
  p<- ggplot(venn_df, aes(B = Positional, C = eQTL, D = Chromatin, A = MAGMA)) +
    geom_venn(fill_color = c ('#F1766D',"#1E90FF","#FF8C00","#4DAF4A"),show_percentage = F,
              set_name_color =c ('#F1766D', "#1E90FF","#FF8C00","#4DAF4A"),
              set_names =c('GWGAS','Positional mapping', 'eQTL mapping','3D Chromatin mapping'),
              set_name_size = 7,text_size = 6) +
    theme_void() +
    scale_x_continuous(limits = c(-3, 3)) +
    scale_y_continuous(limits = c(-3, 3)) +
    coord_fixed()
  p
  
  venn_font <- function(p, font)
  {
    
    grep_grob <- function(gt, lab){
      which(sapply(gt, function(x) grepl(lab, x$name)))
    }
    
    p2 <- ggplot_gtable(ggplot_build(p))
    mygrobs <- p2$grobs
    panel_grob <- mygrobs[[grep_grob(mygrobs, "panel")]]
    venn_grob <- panel_grob$children[[grep_grob(panel_grob$children, "venn")]]
    text_grob <- venn_grob$children[grep_grob(venn_grob$children, "text")]
    text_grob <- do.call(grid::gList, 
                         lapply(text_grob, 
                                function(x) {x$gp$fontfamily <- font; 
                                x}))
    venn_grob$children[grep_grob(venn_grob$children, "text")] <- text_grob
    panel_grob$children[[grep_grob(panel_grob$children, "venn")]] <- venn_grob
    mygrobs[[grep_grob(mygrobs, "panel")]] <- panel_grob
    p2$grobs <- mygrobs
    grid::grid.newpage()
    grid::grid.draw(p2)
    return(p2)
  }
  
  p2 <- venn_font(p, font = "Times New Roman")
  p2
  ggsave(paste0(folder_path, trait, '_vnn.png'), p2, dpi=600,height =10,width = 12 )
  ggsave(paste0(folder_path, trait, '_vnn.pdf'),p2, device = cairo_pdf,height =10,width = 12 )
  
  
  library(openxlsx)
  
  # 创建 workbook
  wb <- createWorkbook()
  
  addWorksheet(wb, "POS")
  addWorksheet(wb, "eQTL")
  addWorksheet(wb, "CHR3D")
  addWorksheet(wb, "MAGMA")
  
  writeData(wb, "POS", pos)
  writeData(wb, "eQTL", eqtl)
  writeData(wb, "CHR3D", chr_3d)
  writeData(wb, "MAGMA", magma)
  
  saveWorkbook(wb,file = paste0(folder_path, trait, '_genes_summary.xlsx'),overwrite = TRUE)
  
}

# CPASSOC
pos <- fread('result/F3_FUMA/CPASSOC/hip_position/hip_position_genes.txt')
eqtl <- fread('result/F3_FUMA/CPASSOC/hip_eQTL/hip_eQTL_genes.txt')
magma <- fread('result/F3_FUMA/CPASSOC/hip_eQTL/hip_magma.genes.out')
magma <- magma %>%
  filter(P <= 0.05 / nrow(magma))

names(pos)
names(eqtl)
names(magma)


library(openxlsx)

# 创建 workbook
wb <- createWorkbook()

addWorksheet(wb, "POS")
addWorksheet(wb, "eQTL")
addWorksheet(wb, "MAGMA")

writeData(wb, "POS", pos)
writeData(wb, "eQTL", eqtl)
writeData(wb, "MAGMA", magma)

saveWorkbook(wb,file = 'result/F3_FUMA/CPASSOC/CPASSOC_genes_summary.xlsx',overwrite = TRUE)


#N====
# 基因 =====
library(survival)
traits <- c('double_right', 'single_right_hip', 'single_right_knee')
for (trait in traits) {
  folder_path <- '/path/to/epi/result/F3_FUMA/'
  
  pos <- fread(paste0(folder_path, trait, '/POS/genes.txt'))
  eqtl <- fread(paste0(folder_path, trait, '/eQTL/genes.txt'))
  magma <- fread(paste0(folder_path, trait, '/POS/magma.genes.out'))
  magma <- magma %>%
    filter(P <= 0.05 / nrow(magma))
  
  all_genes <- Reduce(union, c(pos$symbol, eqtl$symbol, magma$SYMBOL))
  all_ensg <- Reduce(union, c(pos$ensg, eqtl$ensg, magma$GENE))
  genes <- data.frame(
    genes = all_ensg,
    symbol = all_genes
  )
  cat(trait)
  cat("\n")
  cat(all_genes, sep = ",")
  cat("\n")
  write.csv(genes, paste0(folder_path, trait, 'genes.csv'), row.names = F)
}



library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(reactome.db)
library(gground)
library(ggprism)
library(tidyverse)

for (trait in traits) {
  folder_path <- '/path/to/epi/result/F3_FUMA/'
  
  pos <- fread(paste0(folder_path, trait, '/POS/genes.txt'))
  eqtl <- fread(paste0(folder_path, trait, '/eQTL/genes.txt'))
  magma <- fread(paste0(folder_path, trait, '/POS/magma.genes.out'))
  magma <- magma %>%
    filter(P <= 0.05 / nrow(magma))
  
  all_genes <- Reduce(union, c(pos$ensg, eqtl$ensg, magma$ensg))
  
  cat(all_genes, sep = "\n")
  # 去重
  interactors <- gsub("\\..*$", "", all_genes)
  interactors <- unique(interactors)
  
  # 2. ENSEMBL -> ENTREZID
  id <- bitr(interactors, fromType = "ENSEMBL", toType   = "ENTREZID", OrgDb    = org.Hs.eg.db)
  
  # 已经成功匹配的 Ensembl ID
  matched <- unique(id$ENSEMBL)
  
  # 找出没有匹配上的 Ensembl ID
  fail <- setdiff(interactors, matched)
  
  cat("Total Ensembl:", length(interactors), "\n")
  cat("Matched:", length(matched), "\n")
  cat("Failed:", length(fail), "\n")
  
  # 最终用于富集分析的 Entrez ID
  id_final <- id %>% distinct(ENSEMBL, ENTREZID)
  
  # 如果只需要 Entrez ID
  gene_entrez <- unique(id_final$ENTREZID)
  
  # 5. GO enrichment
  go <- enrichGO(
    gene          = gene_entrez,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  
  GO <- go@result
  
  write.csv(GO, paste0("result/F3_FUMA富集分析/GO_", trait, '.csv'), row.names = FALSE)
  
  # 6. KEGG enrichment
  R.utils::setOption("clusterProfiler.download.method", "auto")
  
  kegg <- enrichKEGG(
    gene          = gene_entrez,
    organism      = "hsa",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  KEGG <- kegg@result
  write.csv(KEGG, paste0("result/F3_FUMA富集分析/KEGG_", , trait, '.csv'), row.names = FALSE)
}
# 韦恩图 trait all gene ====

double_genes <- fread('result/F3_FUMA/double_rightgenes.csv')
hip_genes <- fread('result/F3_FUMA/single_right_hipgenes.csv')
knee_genes <- fread('result/F3_FUMA/single_right_kneegenes.csv')


library(ggvenn)
unique_gene <- unique(c(knee_genes$genes, hip_genes$genes, double_genes$genes))
common_genes <- Reduce(intersect, list(knee_genes$genes, hip_genes$genes, double_genes$genes ))
venn_df <- data.frame(gene_ensg = unique_gene,
                      Knee = unique_gene %in% knee_genes$genes,
                      Hip = unique_gene %in% hip_genes$genes,
                      Knee_Hip = unique_gene %in% double_genes$genes)

p <- ggplot(venn_df, aes(A = Knee, B = Hip, C = Knee_Hip)) +
  geom_venn(fill_color = c ('#006affff',"#1dcc2fff","#FD3E5Cff"),show_percentage = F,
            set_name_color =c ('#006affff', "#1dcc2fff","#FD3E5Cff"),
            set_names =c('Knee','Hip', 'Knee - Hip'),
            set_name_size = 7,text_size = 6) +
  theme_void() +
  scale_x_continuous(limits = c(-3, 3)) +
  scale_y_continuous(limits = c(-3, 3)) +
  coord_fixed()
p

library(dplyr)
library(ggvenn)
library(ggplot2)

A <- knee_genes$symbol
B <- hip_genes$symbol
C <- double_genes$symbol

get_genes <- function(region) {
  switch(region,
         "A"   = setdiff(A, union(B, C)),
         "B"   = setdiff(B, union(A, C)),
         "C"   = setdiff(C, union(A, B)),
         "AB"  = setdiff(intersect(A, B), C),
         "AC"  = setdiff(intersect(A, C), B),
         "BC"  = setdiff(intersect(B, C), A),
         "ABC" = Reduce(intersect, list(A, B, C))
  )
}

regions <- c("A", "B", "C", "AB", "AC", "BC", "ABC")

# 区域圆心（按 ggvenn 默认 3 集合布局估算）
center <- data.frame(
  region = c("A", "B", "C", "AB", "AC", "BC", "ABC"),
  x = c(-1.7, 1.7, 0.0, 0.0, -1.2, 1.2, 0.0),
  y = c(1.3, 1.3, -1.3, 1.8, -0.2, -0.2, 0.5)
)

# 每个基因一条记录
gene_label <- do.call(rbind, lapply(regions, function(r) {
  g <- get_genes(r)
  if (length(g) == 0) return(NULL)
  data.frame(region = r, gene = g, stringsAsFactors = FALSE)
})) %>%
  left_join(center, by = "region")

max_per_col <- 4      # 每列最多放几个基因，按需调
col_gap     <- 0.7     # 列间距（坐标单位）
row_step    <- 0.16    # 行距

gene_label <- gene_label %>%
  group_by(region) %>%
  mutate(
    n     = n(),
    # 分列：前 max_per_col 个在第 1 列，接下来在第 2 列……
    col   = ceiling(row_number() / max_per_col),
    ncol  = max(col),
    # 列内行号：每列从 1 重新开始
    row   = row_number() - (col - 1) * max_per_col,
    # 每列的实际长度（最后一列可能不满）
    col_n = pmin(max_per_col, n - (col - 1) * max_per_col),
    # x：以圆心为中心，向左右展开
    x_pos = x + (col - (ncol + 1) / 2) * col_gap,
    # y：每列内部上下居中
    y_pos = y + (col_n - 2 * row + 1) * row_step / 2
  ) %>%
  ungroup()

p <- ggplot(venn_df, aes(A = Knee, B = Hip, C = Knee_Hip)) +
  geom_venn(fill_color = c('#006affff', "#1dcc2fff", "#FD3E5Cff"),
            show_percentage = FALSE,
            set_name_color = c('#006affff', "#1dcc2fff", "#FD3E5Cff"),
            set_names = c('Knee', 'Hip', 'Knee - Hip'),
            set_name_size = 7, text_size = 6) +
  geom_text(data = gene_label,
            aes(x = x_pos, y = y_pos, label = gene),
            inherit.aes = FALSE,
            size = 3.8, family = "Times New Roman", fontface = "bold",
            hjust = 0.5, vjust = 0.5) +
  theme_void() +
  scale_x_continuous(limits = c(-4, 4)) +
  scale_y_continuous(limits = c(-4, 4)) +
  coord_fixed()

p



p2 <- venn_font(p, font = "Times New Roman")
p2
ggsave('result/F3_FUMA/vnn_trait.png', p2, dpi=600,height =10,width = 12 )
ggsave('result/F3_FUMA/vnn_trait.pdf',p2, device = cairo_pdf,height =10,width = 12 )

# N =====
# PRS 疾病死亡====
# 疾病数据
disease_data <- fread("shared_project/data/disease_death/incident_death_disease_v0.csv")

# 协变量
cov_data <- fread("shared_project/data/sex_age_race.csv")
cov_data <- cov_data[,c('eid', 'age_i0', 'education_class')]
# 基因协变量
gene_cov <- fread('data/gwas_cov.csv')
gene_cov <- subset(gene_cov, select = -eid1)
names(gene_cov)[names(gene_cov) == "eid"] <- "IID"

# eid
eid <- fread('shared_project/data/eid_241031.csv')
colnames(eid) <- c('IID','eid')


#  PRS 协变量
base_covariates <- c("age_i0", "sex", "education_class", "batch", "PC1", "PC2", 
                     "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", 
                     "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", 
                     "PC17", "PC18", "PC19", "PC20")

traits <- c("single_right_knee", "single_right_hip", "double_right")

result_list <- list()
for (trait in traits) {
  prs_cs <- fread(paste0("result/F3_PRS/PRS_CS/",trait,"_PRS_non-GWAS.csv"))
  names(prs_cs)[names(prs_cs) == 'PRS'] <- 'PRS_CS'
  
  prs_ct <- fread(paste0("result/F3_PRS/PRS_CT/",trait,"_clump_0.05_PRS.txt"))
  names(prs_ct)[names(prs_ct) == 'PRS'] <- 'PRS_CT'
  

  data <- eid %>%
    inner_join(prs_cs, by = "IID") %>%
    inner_join(prs_ct, by = "IID") %>%
    inner_join(gene_cov, by = "IID")
  
  data <- data %>%
    inner_join(cov_data, by = "eid") %>%
    inner_join(disease_data, by = "eid")
  
  # 疾病列表
  disease_names <- c(
    "death", "CVD", "respiratory", "digestive", "endocrine", "neurological", 
    "psychiatric", "musculoskeleta"
  )
  
  for (exposure in c('PRS_CS', 'PRS_CT')) {
    for (disease in disease_names) {
      
      time_var <- paste0(disease, "_time")
      
      # 检查变量
      if (!disease %in% names(data)) {
        return(NULL)
      }
      if (!time_var %in% names(data)) {
        return(NULL)
      }
      
      # 当前疾病分析数据
      vars_needed <- c(disease, time_var, exposure, base_covariates)
      vars_needed <- unique(vars_needed)
      dat <- data %>%
        filter(complete.cases(across(all_of(vars_needed))))
      
      # 样本量
      N <- nrow(dat)
      events <- sum(dat[[disease]] == 1, na.rm = TRUE)
      if (N == 0) {
        return(NULL)
      }
      
      if (events < 10) {
        warning(disease, ": events < 10")
        return(NULL)
      }
      
      # Cox formula
      covariate_formula <- paste(c(exposure, base_covariates), collapse = " + ")
      formula <- as.formula(paste0("Surv(", time_var, ", ", disease, ") ~ ", covariate_formula))
      
      # Cox model
      fit <- coxph(formula, data = dat)
      
      
      if (is.null(fit)) {
        return(NULL)
      }
      
      # Summary
      fit_summary <- summary(fit)
      coef_table <- fit_summary$coefficients
      conf_table <- fit_summary$conf.int
      
      if (!exposure %in% rownames(coef_table)) {
        return(NULL)
      }
      
      # Age gap result
      result <- data.frame(
        trait = trait,
        exposure = exposure,
        outcome = disease,
        N = N,
        events = events,
        beta = coef_table[exposure, "coef"],
        SE = coef_table[exposure, "se(coef)"],
        HR = conf_table[exposure, "exp(coef)"],
        CI_low = conf_table[exposure, "lower .95"],
        CI_high = conf_table[exposure, "upper .95"],
        P = coef_table[exposure, "Pr(>|z|)"],
        stringsAsFactors = FALSE
      )
      result_list <- rbind(result_list,result)
    }
  }
}

write.csv(result_list, 'result/F3_PRS/all_PRS.csv', row.names = F)

result_list <- fread('result/F3_PRS/all_PRS.csv')

plot_data <- result_list %>%
  filter(grepl('CS', exposure)) %>%
  mutate(
    trait = case_when(
      grepl("^double_right", trait) ~ "Right hip and knee",
      grepl("^single_right_hip", trait) ~ "Right hip",
      grepl("^single_right_knee", trait) ~ "Right knee",
      TRUE ~ trait
    ),
    outcome = factor(
      outcome, 
      levels = c("death", "CVD", "respiratory", "digestive", "endocrine", 
                 "neurological", "psychiatric", "musculoskeleta")
    )
  ) %>%
  mutate(
    trait = factor(trait, levels = c("Right knee", "Right hip", "Right hip and knee")),
    outcome = dplyr::recode(
      outcome,
      death = "All-cause mortality",
      CVD = "Cardiovascular disease",
      respiratory = "Respiratory disease",
      digestive = "Digestive disease",
      endocrine = "Endocrine disease",
      neurological = "Neurological disease",
      psychiatric = "Psychiatric disease",
      musculoskeleta = "Musculoskeletal disease"
    )
  ) %>%
  mutate(
    trait = dplyr::recode(
      trait,
      `Right knee` = "Knee",
      `Right hip` = "Hip",
      `Right hip and knee` = "Knee - Hip"
    )
  )

p <- ggplot(plot_data, aes(x = outcome, y = HR, color = trait)) +
  geom_hline(yintercept = 1,colour = "grey20",linetype = "dashed",linewidth = 0.6) +
  geom_linerange(aes(ymin = CI_low, ymax = CI_high),
                 position = position_dodge(0.6),linewidth = 2) +
  geom_point(position = position_dodge(0.6), shape = 15, size = 5 ) +
  #scale_y_continuous(limits = c(0.99, 1.08), breaks = seq(1, 1.08, 0.02)) +
  scale_color_manual(
    values = c("#006affff","#1dcc2fff","#FD3E5Cff"), name = NULL) +
  labs(x=NULL, y="Hazard ratio (95% CI)", title = NULL) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 3.8))) +
  theme_classic(base_size = 20) +
  theme(
    axis.line = element_line(linewidth = 1),
    axis.ticks = element_line(linewidth = 1),
    axis.ticks.length = unit(0.25, "cm"),
    axis.title = element_text(size=20, family = "Times New Roman", face = "bold"),
    axis.text = element_text(family = "Times New Roman", face = "bold"),
    axis.text.x = element_text(size=20, angle = 60, hjust = 1),
    legend.text = element_text(family = "Times New Roman", face = "bold"),
    legend.position = "top"
  )

p

ggsave('result/F3_PRS/PRS_disease.png', p, width = 9, height = 7,  dpi = 600)




traits <- c('double_right', 'single_right_hip', 'single_right_knee')

for (trait in traits) {
  prs_cs <- fread(paste0("result/F3_PRS/PRS_CS/",trait,"_PRS_GWAS.csv"))
  names(prs_cs)[names(prs_cs) == 'PRS'] <- 'PRS_CS'
  
  prs_ct <- fread(paste0("result/F3_PRS/PRS_CT/",trait,"_clump_0.05_PRS.txt"))
  names(prs_ct)[names(prs_ct) == 'PRS'] <- 'PRS_CT'
  
  pred <- fread(paste0('data/', trait, '_oof.csv'))
  
  data <- eid %>%
    inner_join(prs_cs, by = "IID") %>%
    inner_join(prs_ct, by = "IID") 
  
  data <- data %>%
    inner_join(pred, by = "eid") 
  
  plot(data$PRS_CT, data$age_gap_adjusted)
  plot(data$PRS_CS, data$age_gap_adjusted)
  plot(data$PRS_CS, data$PRS_CT)
}


# PRS biomaker ====

# 指标
cor <- fread('data/cor_data_i0.csv')
names(cor)

# 协变量
cov_data <- fread("shared_project/data/sex_age_race.csv")
cov_data <- cov_data[,c('eid', 'age_i0', 'education_class')]
# 基因协变量
gene_cov <- fread('data/gwas_cov.csv')
gene_cov <- subset(gene_cov, select = -eid1)
names(gene_cov)[names(gene_cov) == "eid"] <- "IID"

# eid
eid <- fread('shared_project/data/eid_241031.csv')
colnames(eid) <- c('IID','eid')


# 协变量
base_covariates <- c("age_i0", "sex", "education_class", "batch", "PC1", "PC2", 
                     "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", 
                     "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", 
                     "PC17", "PC18", "PC19", "PC20")

traits <- c("single_right_knee", "single_right_hip", "double_right")

result_list <- list()
for (trait in traits) {
  prs_cs <- fread(paste0("result/F3_PRS/PRS_CS/",trait,"_PRS_non-GWAS.csv"))
  names(prs_cs)[names(prs_cs) == 'PRS'] <- 'PRS_CS'
  
  data <- eid %>%
    inner_join(prs_cs, by = "IID") %>%
    inner_join(gene_cov, by = "IID")
  
  data <- data %>%
    inner_join(cov_data, by = "eid") %>%
    inner_join(cor, by = "eid")
  
  # biomakers
  markers <- c('father_death_age','mother_death_age','tl','SII',
               'CRP','frailty_score','BMD','grip_min','grip_max')
  
  for (marker in markers) {
    
    exposure <- 'PRS_CS'
    
    # 当前疾病分析数据
    vars_needed <- c(marker, exposure, base_covariates)
    vars_needed <- unique(vars_needed)
    dat <- data %>%
      filter(complete.cases(across(all_of(vars_needed))))
    
    # Linear regression
    covariate_formula <- paste(c(exposure, base_covariates), collapse = " + ")
    formula <- as.formula(paste0(marker, ' ~ ', covariate_formula))
    
    # Linear regression model
    fit <- lm(formula, data = dat)
    
    # Summary
    coef_result <- summary(fit)$coefficients[exposure, ]
    ci <- confint(fit, exposure)
    
    # Age gap result
    result <- data.frame(
      trait = trait,
      exposure = exposure,
      outcome = marker,
      N = nrow(dat),
      beta = coef_result[1],
      SE = coef_result[2],
      CI_low = ci[1],
      CI_high = ci[2],
      P = coef_result[4],
      stringsAsFactors = FALSE
    )
    result_list <- rbind(result_list,result)
  }
}


write.csv(result_list, 'result/F3_PRS/all_PRS_biomaker.csv', row.names = F)

result_list <- fread('result/F3_PRS/all_PRS_biomaker.csv')

plot_data <- result_list %>%
  filter(grepl('CS', exposure)) %>%
  #filter(!grepl('death', outcome))  %>%
  filter(!outcome %in% c('CRP', 'SII', 'grip_min')) %>%
  mutate(
    trait = case_when(
      grepl("^single_right_knee", trait) ~ "Knee",
      grepl("^single_right_hip", trait) ~ "Hip",
      grepl("^double_right", trait) ~ "Knee - Hip",
      TRUE ~ trait
    )
  ) %>%
  mutate(
    trait = factor(trait, levels = c("Knee", "Hip", "Knee - Hip"))
  ) 


# 4. 指标顺序
marker_mapping <- c(
  "father_death_age" = "Father's age at death",
  "mother_death_age" = "Mother's age at death",
  "tl" = "Telomere length",
  "frailty_score" = "Frailty score",
  "BMD" = "Heel bone nineral density",
  "grip_max" = "Hand grip strength"
)

plot_data <- plot_data %>%
  mutate(outcome = factor(outcome, levels = names(marker_mapping), labels = unname(marker_mapping))) %>%
  mutate(outcome = factor(outcome,
                         levels =  c("Father's age at death", "Mother's age at death", 
                                         "Telomere length", 'Heel bone nineral density', 
                                         "Hand grip strength","Frailty score")))



library(ggh4x) 
p <- ggplot(plot_data, aes(x = trait, y = beta, color = trait)) +
  geom_hline(yintercept = 0,colour = "grey20",linetype = "dashed", linewidth = 0.6) +
  geom_linerange(aes(ymin = CI_low, ymax = CI_high), linewidth = 2) +
  geom_point(shape = 15, size = 5) +
  facet_wrap(~ outcome, ncol = 2, strip.position = "top",scales = "free_y") +
  facetted_pos_scales(
    y = list(
      outcome == "Father's age at death" ~ scale_y_continuous(limits = c(-0.8, 0.02), breaks = seq(-0.8, 0,  0.2)),
      outcome == "Mother's age at death" ~ scale_y_continuous(limits = c(-0.65, 0.02), breaks = seq(-0.6, 0,  0.2)),
      outcome == "Telomere length" ~ scale_y_continuous(limits = c(-0.005, 0.001), breaks = seq(-0.005, 0.001,  0.001)),
      outcome == "Heel bone nineral density" ~ scale_y_continuous(limits = c(-0.022, 0.001), breaks = seq(-0.022, 0, 0.011)),
      outcome == "Hand grip strength" ~ scale_y_continuous(limits = c(-0.33, 0.02), breaks = seq(-0.33, 0, 0.11)),
      outcome == "Frailty score" ~ scale_y_continuous(limits = c(-0.002, 0.033), breaks = seq(0, 0.033, 0.011)),
      TRUE ~ scale_y_continuous()
    )
  ) +
  scale_color_manual(values = c("#006affff","#1dcc2fff","#FD3E5Cff")) +
  labs(x = NULL, y = "Regression coefficient (95% CI)") +
  theme_classic(base_size = 20)+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(family = "Times New Roman", face = 'bold', color = 'black'),
    axis.text.y = element_text(size = 14, family = "Times New Roman", face = 'bold', color = 'black'),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    axis.line = element_blank(),
    legend.text = element_text(family = "Times New Roman", face = 'bold'),
    legend.key.size = unit(1, "cm"),
    panel.border = element_rect(linewidth = 1.5, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(size = 16, family = "Times New Roman", face = 'bold', color = 'black'),
    strip.background = element_rect(color = "black", linewidth = 1.2,fill = "grey90"),
    legend.position = 'top'
  )
p
ggsave('result/F3_PRS/PRS_biomaker.pdf', p, width = 8, height = 10,  device = cairo_pdf)
ggsave('result/F3_PRS/PRS_biomaker.png', p, width = 8, height = 10,  dpi = 600)
# PRS 蛋白====

# 所有蛋白数据
protein = fread('/path/to/shared_project/data/res_data/all_protein_scaled.csv')
names(protein)

# 协变量
cov_data <- fread("shared_project/data/sex_age_race.csv")
cov_data <- cov_data[,c('eid', 'age_i0', 'education_class')]
# 基因协变量
gene_cov <- fread('data/gwas_cov.csv')
gene_cov <- subset(gene_cov, select = -eid1)
names(gene_cov)[names(gene_cov) == "eid"] <- "IID"

# eid
eid <- fread('shared_project/data/eid_241031.csv')
colnames(eid) <- c('IID','eid')

# 协变量
base_covariates <- c("age_i0", "sex", "education_class", "batch", "PC1", "PC2", 
                     "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", 
                     "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", 
                     "PC17", "PC18", "PC19", "PC20")

traits <- c("single_right_knee", "single_right_hip", "double_right")

result_list <- list()
for (trait in traits) {
  prs_cs <- fread(paste0("result/F3_PRS/PRS_CS/",trait,"_PRS_non-GWAS.csv"))
  names(prs_cs)[names(prs_cs) == 'PRS'] <- 'PRS_CS'
  
  data <- eid %>%
    inner_join(prs_cs, by = "IID") %>%
    inner_join(gene_cov, by = "IID")
  
  data <- data %>%
    inner_join(cov_data, by = "eid") %>%
    inner_join(protein, by = "eid")
  
  # proteins
  proteins <- names(protein)[names(protein) != "eid"]
  
  for (protein_var in proteins) {
    
    exposure <- 'PRS_CS'
    
    vars_needed <- c(protein_var, exposure, base_covariates)
    vars_needed <- unique(vars_needed)
    dat <- data %>%
      filter(complete.cases(across(all_of(vars_needed))))
    
    # Linear regression
    covariate_formula <- paste(c(exposure, base_covariates), collapse = " + ")
    formula <- as.formula(paste0(protein_var, ' ~ ', covariate_formula))
    
    # Linear regression model
    fit <- lm(formula, data = dat)
    
    # Summary
    coef_result <- summary(fit)$coefficients[exposure, ]
    ci <- confint(fit, exposure)
    
    # Age gap result
    result <- data.frame(
      trait = trait,
      exposure = exposure,
      outcome = protein_var,
      N = nrow(dat),
      beta = coef_result[1],
      SE = coef_result[2],
      CI_low = ci[1],
      CI_high = ci[2],
      P = coef_result[4],
      stringsAsFactors = FALSE
    )
    result_list <- rbind(result_list,result)
  }
}


write.csv(result_list, 'result/F3_PRS/all_PRS_protein.csv', row.names = F)


result_list <- fread('result/F3_PRS/all_PRS_protein.csv')
# 计算 -log10(P)
plot_data <- result_list %>%
  group_by(trait) %>%
  mutate(
    logP = -log10(P),
    significance = case_when(
      P < 0.05 / (nrow(result_list) / 3) ~ trait, 
      TRUE ~ "Not Significant"),) %>%
  ungroup() %>%
  mutate(trait = case_when(
    grepl("single_right_knee", trait) ~ "Knee",
    grepl("single_right_hip", trait)  ~ "Hip",
    grepl("double_right", trait)      ~ "Knee - Hip",
    TRUE ~ NA_character_
  ),
  significant = P < 0.05 / (nrow(result_list) / 3))


label_data <- plot_data %>%
  filter(significant, !is.na(trait)) %>%
  distinct(trait, outcome, .keep_all = TRUE) %>%
  group_by(trait) %>%
  arrange(beta, .by_group = TRUE) %>%
  filter(row_number() == 1 | row_number() == n()) %>%
  ungroup()


library(ggbreak)
library(ggrepel)
p1 <- ggplot(plot_data, aes(x = beta, y = logP)) +
  geom_hline(yintercept = -log10(0.05 / (nrow(result_list) / 3)), linetype = "dashed", linewidth = 0.7, color = "black") +
  # 1. 不显著的点：灰色填充，黑色边框
  geom_point(data = plot_data %>% filter(P >= 0.05 / (nrow(result_list) / 3)), shape = 21, fill = "grey90", 
             color = "grey20", size = 3.5, alpha = 0.7, stroke = 0.1) +
  # 2. 显著的点：按 trait 填充，黑色边框
  geom_point(data = plot_data %>% filter(P < 0.05 / (nrow(result_list) / 3)), aes(fill = trait), shape = 21, 
             color = "black", size = 3.5, stroke = 0.5) +
  scale_fill_manual(values = c("Knee" = "#006affff", "Hip" = "#1dcc2fff", "Knee - Hip" = "#FD3E5Cff")) +
  scale_color_manual(values = c("Knee" = "#006affff", "Hip" = "#1dcc2fff", "Knee - Hip" = "#FD3E5Cff")) +
  geom_text_repel(data = label_data, aes(label = outcome, color = trait),size = 4.5,
                  family = "Times New Roman", fontface = "bold", box.padding   = 0.6, 
                  point.padding = 0.4, min.segment.length = 2,  segment.color = "grey30", 
                  segment.size  = 0.4, max.overlaps  = Inf, show.legend   = F) +
  guides(fill = guide_legend(override.aes = list(size = 6, shape = 21, color = "black"))) +
  labs(x = "Standardized coefficient", y = expression(bold(-log[bold("10")]) ~ bold('(') ~ bolditalic('P') ~ bold('value)'))) +
  theme_classic() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(size = 22, family = "Times New Roman", face='bold', color='black'),  
    axis.title = element_text(size = 22,family = "Times New Roman", face='bold', color='black'),  
    axis.text.x = element_text(size = 18,family = "Times New Roman", face='bold', color='black'),
    axis.text.y = element_text(size = 18,family = "Times New Roman", face='bold', color='black'),
    legend.text = element_text(size = 18, family = "Times New Roman", face = 'bold', color = 'black' ),
    legend.title = element_blank()
  ) 


p1

p1 <- p1 +
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 100, 105), limits = c(0, 105), expand = c(0, 0)) +
  scale_y_break(breaks = c(50, 100), scales = 0.1, space = 0.25) +
  theme(
    axis.line.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank()
  )

p1
ggsave('result/F3_PRS/PRS_protein.png', p1, width = 9, height = 7,  dpi = 600)
ggsave('result/F3_PRS/PRS_protein.pdf', p1, width = 9, height = 7,  device = cairo_pdf)
# PRS 代谢====

meta = fread('/path/to/shared_project/data/res_data/all_meta_scaled.csv')
names(meta)

# 协变量
cov_data <- fread("shared_project/data/sex_age_race.csv")
cov_data <- cov_data[,c('eid', 'age_i0', 'education_class')]
# 基因协变量
gene_cov <- fread('data/gwas_cov.csv')
gene_cov <- subset(gene_cov, select = -eid1)
names(gene_cov)[names(gene_cov) == "eid"] <- "IID"

# eid
eid <- fread('shared_project/data/eid_241031.csv')
colnames(eid) <- c('IID','eid')

# 协变量
base_covariates <- c("age_i0", "sex", "education_class", "batch", "PC1", "PC2", 
                     "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", 
                     "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", 
                     "PC17", "PC18", "PC19", "PC20")

traits <- c("single_right_knee", "single_right_hip", "double_right")

result_list <- list()
for (trait in traits) {
  prs_cs <- fread(paste0("result/F3_PRS/PRS_CS/",trait,"_PRS_non-GWAS.csv"))
  names(prs_cs)[names(prs_cs) == 'PRS'] <- 'PRS_CS'
  
  data <- eid %>%
    inner_join(prs_cs, by = "IID") %>%
    inner_join(gene_cov, by = "IID")
  
  data <- data %>%
    inner_join(cov_data, by = "eid") %>%
    inner_join(meta, by = "eid")
  
  # metas
  metas <- names(meta)[names(meta) != "eid"]
  
  for (meta_var in metas) {
    
    exposure <- 'PRS_CS'
    
    vars_needed <- c(meta_var, exposure, base_covariates)
    vars_needed <- unique(vars_needed)
    dat <- data %>%
      filter(complete.cases(across(all_of(vars_needed))))
    
    # Linear regression
    covariate_formula <- paste(c(exposure, base_covariates), collapse = " + ")
    formula <- as.formula(paste0(meta_var, ' ~ ', covariate_formula))
    
    # Linear regression model
    fit <- lm(formula, data = dat)
    
    # Summary
    coef_result <- summary(fit)$coefficients[exposure, ]
    ci <- confint(fit, exposure)
    
    # Age gap result
    result <- data.frame(
      trait = trait,
      exposure = exposure,
      outcome = meta_var,
      N = nrow(dat),
      beta = coef_result[1],
      SE = coef_result[2],
      CI_low = ci[1],
      CI_high = ci[2],
      P = coef_result[4],
      stringsAsFactors = FALSE
    )
    result_list <- rbind(result_list,result)
  }
}


write.csv(result_list, 'result/F3_PRS/all_PRS_meta.csv', row.names = F)

result_list <- fread('result/F3_PRS/all_PRS_meta.csv')
result_list
meta_name <- fread('shared_project/data/meta_name.csv')
meta_name <- na.omit(meta_name)
meta_name$`Field ID` <- paste0('p',meta_name$UKB.Field.ID)
meta_name <- meta_name[,c('Field ID', 'Metabolite','Shortened name')]
colnames(meta_name) <- c('outcome', 'Metabolite','Shortened name')
head(meta_name)

result_list <- inner_join(meta_name, result_list)
write.csv(result_list, 'result/F3_PRS/all_PRS_meta.csv', row.names = F)
nrow(result_list) / 3

plot_data <- result_list %>%
  group_by(trait) %>%
  mutate(
    logP = -log10(P),
    significance = case_when(
      P < 0.05 / (nrow(result_list) / 3) ~ trait, 
      TRUE ~ "Not Significant"),) %>%
  ungroup() %>%
  mutate(trait = case_when(
    grepl("single_right_knee", trait) ~ "Knee",
    grepl("single_right_hip", trait)  ~ "Hip",
    grepl("double_right", trait)      ~ "Knee - Hip",
    TRUE ~ NA_character_
  ),
  significant = P < 0.05 / (nrow(result_list) / 3))


label_data <- plot_data %>%
  filter(significant, !is.na(trait)) %>%
  distinct(trait, outcome, .keep_all = TRUE) %>%
  group_by(trait) %>%
  arrange(beta, .by_group = TRUE) %>%
  filter(row_number() == 1 | row_number() == n()) %>%
  ungroup()

library(ggbreak)
p2 <- ggplot(plot_data, aes(x = beta, y = logP)) +
  geom_hline(yintercept = -log10(0.05 / (nrow(result_list) / 3)), linetype = "dashed", linewidth = 0.7, color = "black") +
  # 1. 不显著的点：灰色填充，黑色边框
  geom_point(data = plot_data %>% filter(P >= 0.05 / (nrow(result_list) / 3)), shape = 21, fill = "grey90", 
             color = "grey20", size = 3.5, alpha = 0.7, stroke = 0.1) +
  # 2. 显著的点：按 trait 填充，黑色边框
  geom_point(data = plot_data %>% filter(P < 0.05 / (nrow(result_list) / 3)), aes(fill = trait), shape = 21, 
             color = "black", size = 3.5, stroke = 0.5) +
  scale_fill_manual(values = c("Knee" = "#006affff", "Hip" = "#1dcc2fff", "Knee - Hip" = "#FD3E5Cff")) +
  scale_color_manual(values = c("Knee" = "#006affff", "Hip" = "#1dcc2fff", "Knee - Hip" = "#FD3E5Cff")) +
  geom_text_repel(data = label_data, aes(label = `Shortened name`, color = trait),size = 4.5,
                  family = "Times New Roman", fontface = "bold", box.padding   = 0.6, 
                  point.padding = 0.4, min.segment.length = 2,  segment.color = "grey30", 
                  segment.size  = 0.4, max.overlaps  = Inf, show.legend   = F) +
  guides(fill = guide_legend(override.aes = list(size = 6, shape = 21, color = "black"))) +
  labs(x = "Standardized coefficient", y = expression(bold(-log[bold("10")]) ~ bold('(') ~ bolditalic('P') ~ bold('value)'))) +
  theme_classic() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(size = 22, family = "Times New Roman", face='bold', color='black'),  
    axis.title = element_text(size = 22,family = "Times New Roman", face='bold', color='black'),  
    axis.text.x = element_text(size = 18,family = "Times New Roman", face='bold', color='black'),
    axis.text.y = element_text(size = 18,family = "Times New Roman", face='bold', color='black'),
    legend.text = element_text(size = 18, family = "Times New Roman", face = 'bold', color = 'black' ),
    legend.title = element_blank()
  ) 

p2
p2 <- p2 +
  scale_y_continuous(
    breaks = c(0, 2, 4, 6, 8, 14, 16),
    limits = c(0, 16),
    expand = c(0, 0)
  ) +
  scale_y_break(
    breaks = c(6, 8, 9.5, 13),
    scales = 0.1,
    space = 0.25
  ) +
  theme(
    axis.line.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank()
  )

p2
ggsave('result/F3_PRS/PRS_meta.png', p2, width = 9, height = 7,  dpi = 600)
ggsave('result/F3_PRS/PRS_meta.pdf', p2, width = 9, height = 7,  device = cairo_pdf)

# GWAS catelog ====
library(readxl)
double_catalog <- read_excel("result/F3_Catalog/catalog.xlsx", sheet = "Double")
hip_catalog <- read_excel("result/F3_Catalog/catalog.xlsx", sheet = "Hip")
knee_catalog <- read_excel("result/F3_Catalog/catalog.xlsx", sheet = "Knee")
CPASSOC_catalog <- read_excel("result/F3_Catalog/catalog.xlsx", sheet = "CPASSOC")
## Knee ====
# 柱状图
# 最大 logp

knee_catalog$prop <- knee_catalog$n / knee_catalog$N
knee_catalog$logp <- -log10(as.numeric(knee_catalog$`adjusted P`))
knee_catalog


max_logp <- max(knee_catalog$logp, na.rm = TRUE)
# 每个 GeneSet 一个 Y 位置
gene_order <- knee_catalog %>%
  group_by(GeneSet) %>%
  summarise(max_logp_gene = max(logp, na.rm = TRUE), .groups = "drop") %>%
  arrange(max_logp_gene) %>%
  mutate(y_id = row_number())


# 保留所有重复结果
plot_data <- knee_catalog %>%
  left_join(gene_order %>% dplyr::select(GeneSet, y_id), by = "GeneSet") %>%
  mutate(
    # prop 映射到左侧
    prop_x = -prop * max_logp,
    # logp 保持右侧
    logp_x = logp,
    ymin = y_id - 0.5,
    ymax = y_id + 0.5
  )

base_theme <- theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", face = "bold"),
    axis.title.x = element_text(size = 20, color = "black", face = "bold"),
    axis.text.x = element_text(size = 16, color = "black", face = "bold"),
    axis.text.y = element_text(size = 16, color = "black", face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.8, lineend = "square"),
    plot.margin = margin(0, 0, 0, 0)     # 关键：外边距归零
  )
# 左图
y_min <- min(gene_order$y_id) - 0.5
y_max <- max(gene_order$y_id) + 0.5

p_left <- ggplot(plot_data) +
  geom_rect(aes(xmin = prop_x, xmax = 0, ymin = ymin, ymax = ymax),
            fill = "#EEAD1A", color = "black", linewidth = 0.5) +
  scale_y_continuous(breaks = gene_order$y_id, labels = gene_order$GeneSet,
                     limits = c(y_min, y_max), expand = c(0, 0)) +
  scale_x_continuous(
    limits = c(-max_logp-0.2, 0),
    breaks = c(-max_logp, -0.8*max_logp, -0.6*max_logp, -0.4*max_logp, -0.2*max_logp, 0),
    labels = c("1.0", "0.8", "0.6", "0.4", "0.2", "0"),
    expand = c(0, 0)
  ) +
  labs(x = "Proportion", y = NULL) +
  base_theme

# 右图
p_right <- ggplot(plot_data) +
  geom_rect(aes(xmin = 0, xmax = logp_x, ymin = ymin, ymax = ymax),
            fill = "#1AEEAD", color = "black", linewidth = 0.5) +
  scale_y_continuous(breaks = gene_order$y_id, labels = gene_order$GeneSet,
                     limits = c(y_min, y_max), expand = c(0, 0)) +
  scale_x_continuous(
    limits = c(0, max_logp + 0.1),
    breaks = c(0, 2, 4, 6, 8),
    expand = c(0, 0)
  ) +
  labs(x = "-log10 adjusted P value", y = NULL) +
  base_theme +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank()
  ) 

# 拼接：无间隙
p1 <- p_left + p_right
p1

# 二值矩阵热图
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. 从 knee_catalog$genes 提取所有基因
gene_data <- knee_catalog %>%
  dplyr::select(GeneSet, genes) %>%
  filter(!is.na(genes)) %>%
  separate_rows(genes, sep = ",") %>%
  mutate(Gene = trimws(genes)) %>%
  filter(Gene != "") %>%
  distinct(GeneSet, Gene)

# 检查
head(gene_data)

# 2. 所有 Gene，去重复
all_genes <- unique(gene_data$Gene)

length(all_genes)
head(all_genes)


# 3. GeneSet 顺序使用图1的顺序
gene_order2 <- gene_order %>% dplyr::select(GeneSet, y_id)

# 4. 建立完整的 X × Y 格子
heatmap_data <- expand.grid(GeneSet = gene_order2$GeneSet, 
                            Gene = all_genes, stringsAsFactors = FALSE) %>%
  left_join(gene_data %>% mutate(present = 1), by = c("GeneSet", "Gene")) %>%
  mutate(present = ifelse(is.na(present), 0, present)) %>%
  left_join(gene_order2, by = "GeneSet") %>%
  mutate(
    Gene = factor(Gene, levels = all_genes)
  )
heatmap_data$Gene

p2 <- ggplot(heatmap_data, aes(x = Gene, y = y_id)) +
  # 所有格子：只填色，不画边框
  geom_tile(aes(fill = factor(present))) +
  # 只给 present == 1 的格子加边框（用 subset 数据）
  geom_tile(data = subset(heatmap_data, present == 1), aes(x = Gene, y = y_id), 
            fill = NA, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("0" = "white", "1" = "#AD1AEE"), guide = "none") +
  scale_y_continuous(breaks = gene_order$y_id,
                     labels = gene_order$GeneSet,
                     expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = "Gene", y = NULL) +
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", face = "bold"),
    axis.text.x = element_text(size = 14, color = "black", face = "bold",
                               angle = 60,hjust = 1),
    axis.text.y = element_text(size = 14, color = "black", face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.8, lineend = "square"),
    # 显式去掉所有网格线（虽然 classic 默认已无）
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_fixed(ratio = 1)
p2

p2 <- p2 + theme(axis.text.y  = element_blank())
p <- p1 + p2 
ggsave('result/F3_Catalog/knee.pdf', p, height = 8, width = 22, device = cairo_pdf)

## Hip ====
# 柱状图
# 最大 logp
hip_catalog$prop <- hip_catalog$n / hip_catalog$N
hip_catalog$logp <- -log10(as.numeric(hip_catalog$`adjusted P`))
hip_catalog


max_logp <- max(hip_catalog$logp, na.rm = TRUE)

# 每个 GeneSet 一个 Y 位置
gene_order <- hip_catalog %>%
  group_by(GeneSet) %>%
  summarise(max_logp_gene = max(logp, na.rm = TRUE), .groups = "drop") %>%
  arrange(max_logp_gene) %>%
  mutate(y_id = row_number())


# 保留所有重复结果
plot_data <- hip_catalog %>%
  left_join(gene_order %>% dplyr::select(GeneSet, y_id), by = "GeneSet") %>%
  mutate(
    # prop 映射到左侧
    prop_x = -prop * max_logp,
    # logp 保持右侧
    logp_x = logp,
    ymin = y_id - 0.5,
    ymax = y_id + 0.5
  )

base_theme <- theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", face = "bold"),
    axis.title.x = element_text(size = 20, color = "black", face = "bold"),
    axis.text.x = element_text(size = 16, color = "black", face = "bold"),
    axis.text.y = element_text(size = 16, color = "black", face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.8, lineend = "square"),
    plot.margin = margin(0, 0, 0, 0)     # 关键：外边距归零
  )
#1AEEAD
#AD1AEE
#EEAD1A

y_min <- min(gene_order$y_id) - 0.5
y_max <- max(gene_order$y_id) + 0.5
# 左图
p_left <- ggplot(plot_data) +
  geom_rect(aes(xmin = prop_x, xmax = 0, ymin = ymin, ymax = ymax),
            fill = "#EEAD1A", color = "black", linewidth = 0.5) +
  scale_y_continuous(breaks = gene_order$y_id, labels = gene_order$GeneSet,
                     limits = c(y_min, y_max), expand = c(0, 0)) +
  scale_x_continuous(
    limits = c(-0.5*max_logp-0.2, 0),
    breaks = c(-0.5*max_logp, -0.4*max_logp, -0.3*max_logp, -0.2*max_logp,-0.1*max_logp, 0),
    labels = c("0.5", "0.4", "0.3", "0.2", "0.1", "0"),
    expand = c(0, 0)
  ) +
  labs(x = "Proportion", y = NULL) +
  base_theme

# 右图
p_right <- ggplot(plot_data) +
  geom_rect(aes(xmin = 0, xmax = logp_x, ymin = ymin, ymax = ymax),
            fill = "#1AEEAD", color = "black", linewidth = 0.5) +
  scale_y_continuous(breaks = gene_order$y_id, labels = gene_order$GeneSet,
                     limits = c(y_min, y_max), expand = c(0, 0)) +
  scale_x_continuous(
    limits = c(0, max_logp + 0.1),
    breaks = c(0, 1, 2, 3),
    expand = c(0, 0)
  ) +
  labs(x = "-log10 adjusted P value", y = NULL) +
  base_theme +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank()
  ) 

# 拼接：无间隙
p1 <- p_left + p_right
p1

# 二值矩阵热图
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. 从 knee_catalog$genes 提取所有基因
gene_data <- hip_catalog %>%
  dplyr::select(GeneSet, genes) %>%
  filter(!is.na(genes)) %>%
  separate_rows(genes, sep = ",") %>%
  mutate(Gene = trimws(genes)) %>%
  filter(Gene != "") %>%
  distinct(GeneSet, Gene)

# 检查
head(gene_data)

# 2. 所有 Gene，去重复
all_genes <- unique(gene_data$Gene)

length(all_genes)
head(all_genes)

# 3. GeneSet 顺序使用图1的顺序
gene_order2 <- gene_order %>% dplyr::select(GeneSet, y_id)

# 4. 建立完整的 X × Y 格子
heatmap_data <- expand.grid(GeneSet = gene_order2$GeneSet, 
                            Gene = all_genes, stringsAsFactors = FALSE) %>%
  left_join(gene_data %>% mutate(present = 1), by = c("GeneSet", "Gene")) %>%
  mutate(present = ifelse(is.na(present), 0, present)) %>%
  left_join(gene_order2, by = "GeneSet") %>%
  mutate(
    Gene = factor(Gene, levels = all_genes)
  )
heatmap_data$Gene

p2 <- ggplot(heatmap_data, aes(x = Gene, y = y_id)) +
  # 所有格子：只填色，不画边框
  geom_tile(aes(fill = factor(present))) +
  # 只给 present == 1 的格子加边框（用 subset 数据）
  geom_tile(data = subset(heatmap_data, present == 1), aes(x = Gene, y = y_id), 
            fill = NA, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("0" = "white", "1" = "#AD1AEE"), guide = "none") +
  scale_y_continuous(breaks = gene_order$y_id,
                     labels = gene_order$GeneSet,
                     expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = "Gene", y = NULL) +
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", face = "bold"),
    axis.text.x = element_text(size = 14, color = "black", face = "bold",
                               angle = 60,hjust = 1),
    axis.text.y = element_text(size = 14, color = "black", face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.8, lineend = "square"),
    # 显式去掉所有网格线（虽然 classic 默认已无）
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_fixed(ratio = 1)
p2

p2 <- p2 + theme(axis.text.y  = element_blank())
p <- p1 + p2 
ggsave('result/F3_Catalog/hip.pdf', p, height = 8, width = 22, device = cairo_pdf)


## Knee - Hip ====
# 柱状图
# 最大 logp
double_catalog$prop <- double_catalog$n / double_catalog$N
double_catalog$logp <- -log10(as.numeric(double_catalog$`adjusted P`))
double_catalog

max_logp <- max(double_catalog$logp, na.rm = TRUE)

# 每个 GeneSet 一个 Y 位置
gene_order <- double_catalog %>%
  group_by(GeneSet) %>%
  summarise(max_logp_gene = max(logp, na.rm = TRUE), .groups = "drop") %>%
  arrange(max_logp_gene) %>%
  mutate(y_id = row_number())


# 保留所有重复结果
plot_data <- double_catalog %>%
  left_join(gene_order %>% dplyr::select(GeneSet, y_id), by = "GeneSet") %>%
  mutate(
    # prop 映射到左侧
    prop_x = -prop * max_logp,
    # logp 保持右侧
    logp_x = logp,
    ymin = y_id - 0.5,
    ymax = y_id + 0.5
  )

base_theme <- theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", face = "bold"),
    axis.title.x = element_text(size = 20, color = "black", face = "bold"),
    axis.text.x = element_text(size = 16, color = "black", face = "bold"),
    axis.text.y = element_text(size = 16, color = "black", face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.8, lineend = "square"),
    plot.margin = margin(0, 0, 0, 0)     # 关键：外边距归零
  )

y_min <- min(gene_order$y_id) - 0.5
y_max <- max(gene_order$y_id) + 0.5

# 左图
p_left <- ggplot(plot_data) +
  geom_rect(aes(xmin = prop_x, xmax = 0, ymin = ymin, ymax = ymax),
            fill = "#EEAD1A", color = "black", linewidth = 0.5) +
  scale_y_continuous(breaks = gene_order$y_id, labels = gene_order$GeneSet,
                     limits = c(y_min, y_max), expand = c(0, 0)) +
  scale_x_continuous(
    limits = c(-0.4*max_logp-0.2, 0),
    breaks = c(-0.4*max_logp, -0.3*max_logp, -0.2*max_logp,-0.1*max_logp, 0),
    labels = c("0.4", "0.3", "0.2", "0.1", "0"),
    expand = c(0, 0)
  ) +
  labs(x = "Proportion", y = NULL) +
  base_theme

# 右图
p_right <- ggplot(plot_data) +
  geom_rect(aes(xmin = 0, xmax = logp_x, ymin = ymin, ymax = ymax),
            fill = "#1AEEAD", color = "black", linewidth = 0.5) +
  scale_y_continuous(breaks = gene_order$y_id, labels = gene_order$GeneSet,
                     limits = c(y_min, y_max), expand = c(0, 0)) +
  scale_x_continuous(
    limits = c(0, max_logp + 0.1),
    breaks = c(0, 1,2,3),
    expand = c(0, 0)
  ) +
  labs(x = "-log10 adjusted P value", y = NULL) +
  base_theme +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank()
  ) 

# 拼接：无间隙
p1 <- p_left + p_right
p1

# 二值矩阵热图
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. 从 knee_catalog$genes 提取所有基因
gene_data <- double_catalog %>%
  dplyr::select(GeneSet, genes) %>%
  filter(!is.na(genes)) %>%
  separate_rows(genes, sep = ",") %>%
  mutate(Gene = trimws(genes)) %>%
  filter(Gene != "") %>%
  distinct(GeneSet, Gene)

# 检查
head(gene_data)

# 2. 所有 Gene，去重复
all_genes <- unique(gene_data$Gene)

length(all_genes)
head(all_genes)

# 3. GeneSet 顺序使用图1的顺序
gene_order2 <- gene_order %>% dplyr::select(GeneSet, y_id)

# 4. 建立完整的 X × Y 格子
heatmap_data <- expand.grid(GeneSet = gene_order2$GeneSet, 
                            Gene = all_genes, stringsAsFactors = FALSE) %>%
  left_join(gene_data %>% mutate(present = 1), by = c("GeneSet", "Gene")) %>%
  mutate(present = ifelse(is.na(present), 0, present)) %>%
  left_join(gene_order2, by = "GeneSet") %>%
  mutate(
    Gene = factor(Gene, levels = all_genes)
  )
heatmap_data$Gene

p2 <- ggplot(heatmap_data, aes(x = Gene, y = y_id)) +
  # 所有格子：只填色，不画边框
  geom_tile(aes(fill = factor(present))) +
  # 只给 present == 1 的格子加边框（用 subset 数据）
  geom_tile(data = subset(heatmap_data, present == 1), aes(x = Gene, y = y_id), 
            fill = NA, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("0" = "white", "1" = "#AD1AEE"), guide = "none") +
  scale_y_continuous(breaks = gene_order$y_id,
                     labels = gene_order$GeneSet,
                     expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = "Gene", y = NULL) +
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", face = "bold"),
    axis.text.x = element_text(size = 14, color = "black", face = "bold",
                               angle = 60,hjust = 1),
    axis.text.y = element_text(size = 14, color = "black", face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.8, lineend = "square"),
    # 显式去掉所有网格线（虽然 classic 默认已无）
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_fixed(ratio = 1)
p2

p2 <- p2 + theme(axis.text.y  = element_blank())
p <- p1 + p2 
ggsave('result/F3_Catalog/knee-hip.pdf', p, height = 8, width = 22, device = cairo_pdf)
## CPASSOC ====
# 柱状图
# 最大 logp

CPASSOC_catalog$prop <- CPASSOC_catalog$n / CPASSOC_catalog$N
CPASSOC_catalog$logp <- -log10(as.numeric(CPASSOC_catalog$`adjusted P`))
CPASSOC_catalog


max_logp <- max(CPASSOC_catalog$logp, na.rm = TRUE)
# 每个 GeneSet 一个 Y 位置
gene_order <- CPASSOC_catalog %>%
  group_by(GeneSet) %>%
  summarise(max_logp_gene = max(logp, na.rm = TRUE), .groups = "drop") %>%
  arrange(max_logp_gene) %>%
  mutate(y_id = row_number())


# 保留所有重复结果
plot_data <- CPASSOC_catalog %>%
  left_join(gene_order %>% dplyr::select(GeneSet, y_id), by = "GeneSet") %>%
  mutate(
    # prop 映射到左侧
    prop_x = -prop * max_logp,
    # logp 保持右侧
    logp_x = logp,
    ymin = y_id - 0.5,
    ymax = y_id + 0.5
  )

base_theme <- theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", face = "bold"),
    axis.title.x = element_text(size = 24, color = "black", face = "bold"),
    axis.text.x = element_text(size = 20, color = "black", face = "bold"),
    axis.text.y = element_text(size = 22, color = "black", face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.8, lineend = "square"),
    plot.margin = margin(0, 0, 0, 0)     # 关键：外边距归零
  )
# 左图
y_min <- min(gene_order$y_id) - 0.5
y_max <- max(gene_order$y_id) + 0.5

p_left <- ggplot(plot_data) +
  geom_rect(aes(xmin = prop_x, xmax = 0, ymin = ymin, ymax = ymax),
            fill = "#EEAD1A", color = "black", linewidth = 0.5) +
  scale_y_continuous(breaks = gene_order$y_id, labels = gene_order$GeneSet,
                     limits = c(y_min, y_max), expand = c(0, 0)) +
  scale_x_continuous(
    limits = c(-max_logp-0.2, 0),
    breaks = c(-max_logp, -0.8*max_logp, -0.6*max_logp, -0.4*max_logp, -0.2*max_logp, 0),
    labels = c("1.0", "0.8", "0.6", "0.4", "0.2", "0"),
    expand = c(0, 0)
  ) +
  labs(x = "Proportion", y = NULL) +
  base_theme

# 右图
p_right <- ggplot(plot_data) +
  geom_rect(aes(xmin = 0, xmax = logp_x, ymin = ymin, ymax = ymax),
            fill = "#1AEEAD", color = "black", linewidth = 0.5) +
  scale_y_continuous(breaks = gene_order$y_id, labels = gene_order$GeneSet,
                     limits = c(y_min, y_max), expand = c(0, 0)) +
  scale_x_continuous(
    limits = c(0, max_logp + 0.1),
    breaks = c(0, 2, 4, 6, 8),
    expand = c(0, 0)
  ) +
  labs(x = "-log10 adjusted P value", y = NULL) +
  base_theme +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank()
  ) 


# 二值矩阵热图
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. 从 CPASSOC_catalog$genes 提取所有基因
gene_data <- CPASSOC_catalog %>%
  dplyr::select(GeneSet, genes) %>%
  filter(!is.na(genes)) %>%
  separate_rows(genes, sep = ",") %>%
  mutate(Gene = trimws(genes)) %>%
  filter(Gene != "") %>%
  distinct(GeneSet, Gene)

# 检查
head(gene_data)

# 2. 所有 Gene，去重复
all_genes <- unique(gene_data$Gene)

length(all_genes)
head(all_genes)


# 3. GeneSet 顺序使用图1的顺序
gene_order2 <- gene_order %>% dplyr::select(GeneSet, y_id)

# 4. 建立完整的 X × Y 格子
heatmap_data <- expand.grid(GeneSet = gene_order2$GeneSet, 
                            Gene = all_genes, stringsAsFactors = FALSE) %>%
  left_join(gene_data %>% mutate(present = 1), by = c("GeneSet", "Gene")) %>%
  mutate(present = ifelse(is.na(present), 0, present)) %>%
  left_join(gene_order2, by = "GeneSet") %>%
  mutate(
    Gene = factor(Gene, levels = all_genes)
  )
heatmap_data$Gene

p2 <- ggplot(heatmap_data, aes(x = Gene, y = y_id)) +
  # 所有格子：只填色，不画边框
  geom_tile(aes(fill = factor(present))) +
  # 只给 present == 1 的格子加边框（用 subset 数据）
  geom_tile(data = subset(heatmap_data, present == 1), aes(x = Gene, y = y_id), 
            fill = NA, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("0" = "white", "1" = "#AD1AEE"), guide = "none") +
  scale_y_continuous(breaks = gene_order$y_id,
                     labels = gene_order$GeneSet,
                     expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = "Gene", y = NULL) +
  theme_classic() +
  theme(
    text = element_text(family = "Times New Roman", face = "bold"),
    axis.text.x = element_text(size = 20, color = "black", face = "bold",
                               angle = 60,hjust = 1),
    axis.text.y = element_text(size = 20, color = "black", face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.8, lineend = "square"),
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_fixed(ratio = 1)
p2


# 1. 隐藏 p_right 的 Y 轴标签
p_right <- p_right +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.y  = element_blank()
  )

# 2. 隐藏 p2 的 Y 轴标签
p2 <- p2 + theme(axis.text.y  = element_blank())

# 3. 所有图的外边距归零
p_left <- p_left +
  theme(plot.margin = margin(0, 0, 0, 0))

p_right <- p_right +
  theme(plot.margin = margin(0, 0, 0, 0))

p2 <- p2 +
  theme(plot.margin = margin(0, 0, 0, 0))

# 4. 拼接 p1：两个直方图无间距
p1 <- p_left + p_right +
  plot_layout(widths = c(2, 1))

# 5. 拼接 p1 和 p2：无间距
p <- wrap_plots(
  p1, p2,
  nrow = 1,
  ncol = 2,
  widths = c(1, 1.5)
)

ggsave('result/F3_Catalog/CPASSOC.pdf', p, height = 10, width = 34.8, device = cairo_pdf)
# 各个trait间的遗传相关性====
ld_reg <- fread('result/F3_FUMA/all.rg')
ld_reg
# 定义 trait 映射与顺序
trait_mapping <- c(
  "single_right_knee"     = "Knee",
  "single_right_hip" = "Hip",
  "double_right" = "Knee - Hip"
)
ordered_names <- unname(trait_mapping)  # 映射后的 trait 顺序
# 获取所有 trait
traits <- sort(unique(c(ld_reg$p1, ld_reg$p2)))
# 创建空矩阵并填充 rg
rg_mat <- matrix(NA, nrow = length(traits), ncol = length(traits),
                 dimnames = list(traits, traits))
for (i in 1:nrow(ld_reg)) {
  t1 <- ld_reg$p1[i]
  t2 <- ld_reg$p2[i]
  rg_mat[t1, t2] <- ld_reg$rg[i]
  rg_mat[t2, t1] <- ld_reg$rg[i]  # 双向填充
}
# 将行列名改为映射后的全称
rownames(rg_mat) <- trait_mapping[rownames(rg_mat)]
colnames(rg_mat) <- trait_mapping[colnames(rg_mat)]
# 按映射顺序重排
rg_mat <- rg_mat[ordered_names, ordered_names]
# 清空对角线和上三角
rg_mat[upper.tri(rg_mat, diag = F)] <- NA
# 转成长格式，并过滤 NA
data_long <- reshape2::melt(rg_mat, na.rm = FALSE) %>%
  filter(!is.na(value))  # 只保留有值的格子

# 配色
my_colors <- c("#ebf0b5", "#a0d292", "#44a05c" )
# 绘图
ggplot(data_long, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = round(value, 3)), size = 9, lineheight = 0.9, 
            family = "Times New Roman", colour = "black", fontface = "bold") +
  scale_fill_gradientn(
    colors = my_colors,
    limits = c(0, 1),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
    na.value = "white",
    name = "Correlation coefficients (r)",
    guide = guide_colorbar(
      title.position = "left",        # 标题放到色条左侧
      title.theme = element_text(
        angle = 90,                   # 竖向
        hjust = 0.5,                  # 在色条长度方向上居中
        vjust = 0.5,
        family = "Times New Roman",
        size = 30, face = "bold", colour = "black"
      ),
      barwidth  = unit(0.8, "cm"),
      barheight = unit(13, "cm")      # ← 色条加长
    )
  ) +
  scale_y_discrete(limits = rev(levels(data_long$Var1))) + # 翻转 y 轴
  theme_minimal() +
  theme(
    legend.title = element_text(size = 32, family = "Times New Roman", 
                                colour = "black", face = "bold"),
    legend.text = element_text(family = "Times New Roman", size = 28, face = "bold"),
    axis.text.y = element_text(size = 32, family = "Times New Roman", colour = "black", face = "bold"),
    axis.text.x = element_text(size = 32, family = "Times New Roman", colour = "black", face = "bold"),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white", color = "white"),
    legend.key.size = unit(0.8, "cm")
  ) +
  coord_fixed()

ggsave('result/F3_FUMA/LDSC.pdf', height = 7, width = 10, device = cairo_pdf)
ggsave('result/F3_FUMA/LDSC.png', height = 7, width = 10, dpi = 600)


# MiXeR韦恩图 ====
# Hip--Knee
knee_specific <- 0.1
knee_specific_sd <- 0.1

shared <- 1.2
shared_sd <- 0.2

hip_specific <- 1.0
hip_specific_sd <- 0.3

venn_data <- data_frame(
  data = c('1', '2', '3'),
  Knee = c(T,T,F),
  Hip = c(F,T,T)
)


# ggvenn
p <- ggplot(venn_data, aes(A = Knee, B = Hip)) +
  geom_venn(fill_color = c ('#006affff',"#1dcc2fff"),show_percentage = F,
            set_name_color =c ('#006affff', "#1dcc2fff"),
            set_names =c('Knee','Hip'),
            set_name_size = 7,text_size = 6) +
  theme_void() +
  scale_x_continuous(limits = c(-3, 3)) +
  scale_y_continuous(limits = c(-3, 3)) +
  coord_fixed() 
p
p2 <- venn_font(p, font = "Times New Roman")
p2
ggsave('result/F3_FUMA/vnn_hip_knee.pdf',p2, device = cairo_pdf,height =10,width = 12 )


# Double -- Hip

venn_data <- data_frame(
  data = c('1', '2', '3'),
  Double = c(T,T,F),
  Hip = c(F,T,T)
)

# ggvenn
p <- ggplot(venn_data, aes(A = Double, B = Hip)) +
  geom_venn(fill_color = c ('#FD3E5Cff',"#1dcc2fff"),show_percentage = F,
            set_name_color =c ('#FD3E5Cff', "#1dcc2fff"),
            set_names =c('Knee - Hip','Hip'),
            set_name_size = 7,text_size = 6) +
  theme_void() +
  scale_x_continuous(limits = c(-3, 3)) +
  scale_y_continuous(limits = c(-3, 3)) +
  coord_fixed() 
p
p2 <- venn_font(p, font = "Times New Roman")
p2
ggsave('result/F3_FUMA/vnn_double_hip.pdf',p2, device = cairo_pdf,height =10,width = 12 )

# Double -- Hip
venn_data <- data_frame(
  data = c('1', '2', '3'),
  Double = c(T,T,F),
  Knee = c(F,T,T)
)

# ggvenn
p <- ggplot(venn_data, aes(A = Double, B = Knee)) +
  geom_venn(fill_color = c ('#FD3E5Cff',"#006affff"),show_percentage = F,
            set_name_color =c ('#FD3E5Cff', "#006affff"),
            set_names =c('Knee - Hip','Knee'),
            set_name_size = 7,text_size = 6) +
  theme_void() +
  scale_x_continuous(limits = c(-3, 3)) +
  scale_y_continuous(limits = c(-3, 3)) +
  coord_fixed() 
p
p2 <- venn_font(p, font = "Times New Roman")
p2
ggsave('result/F3_FUMA/vnn_double_knee.pdf',p2, device = cairo_pdf,height =10,width = 12 )


# LeadSNP GWAS ====
traits <- c('double_right', 'single_right_hip', 'single_right_knee')
for (trait in traits) {
  folder_path <- '/path/to/epi/result/F3_FUMA/'
  lead_snp <- fread(paste0(folder_path, trait, '/POS/leadSNPs.txt'))
  print(lead_snp)
  names(lead_snp) <- c("No", "GenomicLocus", "uniqID", "SNP", "CHR","POS", "P","nIndSigSNPs","IndSigSNPs")
  gwas <- fread(paste0('data/GWAS_scPagwas/', trait, '.txt'))
  gwas_lead <- inner_join(gwas, lead_snp)
  
  gwas_lead_out <- gwas_lead %>%
    mutate(
      # Minor allele frequency
      `Minor allele frequency` = pmin(AF1, 1 - AF1),
      # Approximate minor allele count
      `Minor allele count` = round(2 * N * `Minor allele frequency`)
    ) %>%
    dplyr::rename(
      `Lead SNP` = SNP,
      Chromosome = CHR,
      Position = POS,
      `Effect allele` = A1,
      `Other allele` = A2,
      `Effect estimate` = BETA,
      `Standard error` = SE,
      `P value` = P,
      `Number of independent significant SNPs` = nIndSigSNPs,
      `Independent significant SNPs` = IndSigSNPs
    ) %>%
    dplyr::select(
      GenomicLocus,
      `Lead SNP`,
      Chromosome,
      Position,
      `Effect allele`,
      `Other allele`,
      `Minor allele frequency`,
      `Minor allele count`,
      `Effect estimate`,
      `Standard error`,
      `P value`,
      `Number of independent significant SNPs`,
      `Independent significant SNPs`
    )
  
  write.csv(gwas_lead_out, paste0(folder_path, trait, '_leadSNP.csv'), 
            row.names = F, na = "")
}


# MAGMA
traits <- c('double_right', 'single_right_hip', 'single_right_knee')
all_magma <- data.frame()
for (trait in traits) {
  folder_path <- '/path/to/epi/result/F3_FUMA/'
  magma <- data.frame(fread(paste0(folder_path, trait, '/POS/magma.genes.out')))
  magma <- magma %>%
    filter(P <= 0.05 / nrow(magma))
  magma['trait'] <- trait
  all_magma <- rbind(all_magma, magma)
}
write.csv(all_magma, paste0(folder_path, 'all_magma.csv'),  row.names = F)





# N-GWAMA
lead_snp <- fread('result/F3_FUMA/N-GWAMA/POS/leadSNPs.txt')
print(lead_snp)
names(lead_snp) <- c("No", "GenomicLocus", "uniqID", "SNP", "CHR","POS", "P","nIndSigSNPs","IndSigSNPs")
gwas <- fread('data/GWAS_scPagwas/knee_hip.N_weighted_GWAMA.results.txt.gz')
names(gwas)
gwas <- gwas[gwas$SNPID %in% lead_snp$SNP,]
names(gwas)[2] <- 'SNP'


gwas_lead <- inner_join(lead_snp, gwas)
names(gwas_lead)

gwas_lead_out <- gwas_lead %>%
  mutate(
    # Minor allele frequency
    `Minor allele frequency` = pmin(EAF, 1 - EAF)
    # Approximate minor allele count
    #`Minor allele count` = round(2 * N * `Minor allele frequency`)
  ) %>%
  dplyr::rename(
    `Lead SNP` = SNP,
    Chromosome = CHR,
    Position = POS,
    `Effect allele` = EA,
    `Other allele` = OA,
    `Effect estimate` = BETA,
    `Standard error` = SE,
    `P value` = PVAL,
    `Number of independent significant SNPs` = nIndSigSNPs,
    `Independent significant SNPs` = IndSigSNPs
  ) %>%
  dplyr::select(
    GenomicLocus,
    `Lead SNP`,
    Chromosome,
    Position,
    `Effect allele`,
    `Other allele`,
    `Minor allele frequency`,
    #`Minor allele count`,
    `Effect estimate`,
    `Standard error`,
    `P value`,
    `Number of independent significant SNPs`,
    `Independent significant SNPs`
  )

write.csv(gwas_lead_out, 'result/F3_FUMA/N-GWAMA_leadSNP.csv', row.names = F)







# QQ图====
traits <- c('double_right', 'single_right_hip', 'single_right_knee')
for (trait in traits) {
  gwas <- fread(paste0('data/GWAS_scPagwas/', trait, '.txt'))
  d1 = gwas %>% dplyr:: select(Marker = SNP, Chr = CHR, Pos = POS, p=P)
  d1 = d1 %>% drop_na(p)
  nrow(d1)
  library(CMplot)
  # QQ图
  CMplot(d1, plot.type = "q", threshold = 0.05,
         signal.cex = c(0.3,0.3),dpi=300,
         file="pdf",file.name=paste0(trait))
}

# SNP柱状图====
traits <- c('double_right', 'single_right_hip', 'single_right_knee')
for (trait in traits) {
  folder_path <- '/path/to/epi/result/F3_FUMA/'
  
  data <- fread(paste0(folder_path, trait, '/POS/annov.stats.txt'))
  data <- data %>%
    mutate(annot = reorder(annot, -prop),
           label = case_when(
             fisher.P < 0.05/11 ~ "**",  # p < 0.05/11 标注 **
             fisher.P < 0.05 ~ "*",         # p < 0.05 标注 *
             TRUE ~ ""              # 不满足条件时不标注
           ))
  data$annot <- gsub('_',' ', data$annot)
  
  data <- data %>%
    mutate(annot = fct_reorder(factor(annot), -prop))
  # 绘制柱状图
  ggplot(data, aes(x = annot, y = prop, fill = enrichment)) + 
    geom_bar(stat = "identity", width = 0.7) +  # 使用prop确定柱子高度
    geom_text(aes(label = label), vjust = -0.01, size = 5, color = "black") +  # 标注 * 或 **
    scale_fill_gradient2(mid = "white", high = "red", midpoint = 0) +  
    labs(y = 'Proportion', x = '', fill = 'Enrichment value') +  # 修改图例标题
    theme_classic() + 
    theme(
      panel.grid = element_blank(),
      legend.title = element_text(size = 18, family = "Times New Roman", colour = 'black', face = "bold"),
      legend.text = element_text(family = 'Times New Roman', size = 14, face = 'bold'),
      axis.title = element_text(size = 22, family = "Times New Roman", face = 'bold', color = 'black'),  
      axis.text = element_text(size = 18, family = "Times New Roman", face = 'bold', color = 'black'),
      axis.text.x = element_text(size = 18, family = "Times New Roman", face = 'bold', color = 'black', angle = 45, hjust = 1)
    )
  ggsave(paste0(folder_path, trait,'_snp.png'), height=6, width=12,dpi=600)
}

# A ====
traits <- c('double_right', 'single_right_hip', 'single_right_knee')
all_char <- data.frame()

for (trait in traits) {
  ## Summary per genomic risk locus
  chr_data <- fread(paste0('result/F3_FUMA/', trait, '/POS/GenomicRiskLoci.txt'))
  chr_data$genomic_loci <- paste0(chr_data$chr,':', chr_data$start, '-', chr_data$end)
  chr_data <- chr_data[,c('chr', 'start', 'end', 'genomic_loci')]
  chr_data$size <- (chr_data$end - chr_data$start) / 1000
  
  ## SNPs 部分 
  snps <- fread(paste0('result/F3_FUMA/', trait, '/POS/snps.txt'))
  snp_data <- snps %>%
    transmute(snp_id = rsID, chr = chr, pos = pos)
  
  # 区间匹配：SNP 单点落在 chr_data 区间内（within）
  overlap <- snp_data %>%
    inner_join(
      chr_data %>% dplyr::select(chr, start, end, genomic_loci),
      by = "chr",
      relationship = "many-to-many"
    ) %>%
    filter(pos >= start, pos <= end)
  
  
  # 统计每个区间的 SNP 数量
  snp_counts <- overlap %>%
    count(chr, start, end, name = "n_snps")
  
  chr_data <- chr_data %>%
    left_join(snp_counts, by = c("chr", "start", "end")) %>%
    mutate(n_snps = replace_na(n_snps, 0L))
  
  ## mapped genes 
  src_paths <- c(
    pos  = paste0('result/F3_FUMA/', trait, '/POS/genes.txt'),
    eqtl = paste0('result/F3_FUMA/', trait, '/eQTL/genes.txt')
  )
  
  # 不同来源用不同窗口
  WIN_map <- c(pos = 10000L, eqtl = 1000000L)
  
  for (nm in names(src_paths)) {
    WIN <- WIN_map[[nm]]
    
    genes <- fread(src_paths[nm]) %>%
      dplyr::select(ensg, chr, start, end) %>%
      dplyr::rename(g_start = start, g_end = end)
    
    cnt <- genes %>%
      inner_join(
        chr_data %>% dplyr::select(chr, start, end),
        by = "chr",
        relationship = "many-to-many"
      ) %>%
      filter(g_start <= end   + WIN,
             g_end   >= start - WIN) %>%
      count(chr, start, end, name = paste0(nm, "_genes"))
    
    chr_data <- chr_data %>%
      left_join(cnt, by = c("chr", "start", "end"))
  }
  
  gene_cols <- paste0(names(src_paths), "_genes")
  chr_data <- chr_data %>%
    mutate(across(all_of(gene_cols), ~ replace_na(.x, 0L)))
  # genes physically located in loci
  # 手动数#genes physically located in loci柱子高度
  if (grepl("knee", trait)) {
    n_gene_loci <- c(0,3,2,2,3)
  } else if (grepl("hip", trait)) {
    n_gene_loci <- c(0,2,1,0,1,1)
  } else {
    n_gene_loci <- c(2,1,2,1,1)
  }
  
  gene_loci <- data.frame(
    genomic_loci = chr_data$genomic_loci,
    n_gene_loci = n_gene_loci
  )
  
  chr_data <- inner_join(chr_data, gene_loci)
  chr_data$trait <- trait
  all_char <- rbind(all_char, chr_data)
}

write.csv(all_char,'result/F3_FUMA/genomic_loci.csv', row.names = F)

## 绘图====
all_char <- data.table::fread('result/F3_FUMA/genomic_loci.csv')
traits <- c('double_right', 'single_right_hip', 'single_right_knee')
for (trait_val in traits) {
  chr_data <- all_char %>%
    filter(trait == trait_val) %>%
    mutate(
      chr_num = as.integer(chr)            # 若 chr 是 "1","2",…，可直接转整数
    ) %>%
    arrange(chr_num, start) %>%            # 先按染色体，再按坐标
    mutate(
      genomic_loci = factor(genomic_loci, levels = rev(genomic_loci))  # 锁定排序后的顺序
    )
  
  
  plot_data <- tidyr::pivot_longer(chr_data, 
                                   cols = c(n_snps, size,  pos_genes, eqtl_genes, n_gene_loci), 
                                   names_to = "category", 
                                   values_to = "value")
  library(dplyr)
  plot_data <- plot_data %>%
    mutate(category = case_when(
      category == "n_snps" ~ "SNPs",
      category == "size" ~ "Size (kb)",
      category == "pos_genes" ~ "Positional mapping",
      category == "eqtl_genes" ~ "eQTL mapping",
      category == "n_gene_loci" ~ "Genes physically located in loci"
    ))
  plot_data$category <- factor(plot_data$category, 
                               levels = c("Size (kb)", "SNPs", "Genes physically located in loci",
                                          "Positional mapping", "eQTL mapping"))
  # 定义颜色映射
  life_color <- c("#F1766D","#746FB1","#87bade","#96D274", "#FCBB44")
  unique(plot_data$category)
  
  
  ggplot(plot_data, aes(x = genomic_loci, y = value, fill = category)) +
    geom_bar(stat = "identity", position = "dodge") +  
    facet_wrap(~ category, scales = "free_x",nrow = 1) +  # 允许y轴（原x轴）独立
    labs(x = "Genomic Loci", y = "", fill = "Category") +
    scale_fill_manual(values = life_color) +
    theme_bw() + 
    theme(
      legend.position = 'none',
      panel.grid = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      strip.text = element_text(size = 14,family = "Times New Roman", face='bold', color='black'),
      strip.background = element_rect(color = "black", size = 1),  # 加粗分面面板的框线
      panel.border = element_rect(size = 1, color = "black"),
      axis.title = element_text(size = 18,family = "Times New Roman", face='bold', color='black'),  
      axis.text = element_text(size = 12,family = "Times New Roman", face='bold', color='black')
      #axis.text.y = element_text(size = 16,family = "Times New Roman", face='bold', color='black')
      
    ) +
    coord_flip()
  if (grepl("hip", trait_val)) {
    ggsave(paste0('result/F3_FUMA/',trait_val, '_snp_gene.png'), height=3, width=19,dpi = 600)
  } else {
    ggsave(paste0('result/F3_FUMA/',trait_val, '_snp_gene.png'), height=2.5, width=19,dpi = 600)
  }
}


