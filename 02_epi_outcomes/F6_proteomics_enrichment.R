# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 02_epi_outcomes
# 蛋白/代谢组关联与功能富集分析
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
# 蛋白====

# 所有蛋白数据
protein = fread('/path/to/shared_project/data/res_data/all_protein_scaled.csv')
names(protein)

# 随访时间 V0-V2
p53 <- fread('shared_project/data/p53.csv')
p53 <- p53 %>%
  mutate(across(where(is.character), ~ if_else(str_detect(., "^\\s*$"), NA_character_, .))) %>%
  mutate(followup_time = as.numeric(as.Date(p53_i2) - as.Date(p53_i0)))

# V0 age
age_v0 <- fread('shared_project/data/sex_age_race.csv')
age_v0 <- age_v0[,c(1,3)]
# 预测年龄
csv_files <- list.files(path = 'data', 
                        pattern = ".*oof.*\\.csv$", full.names = TRUE)

all_results <- list()

for (csv_file in csv_files) {
  pred_age <- fread(csv_file)
  
  data <- pred_age %>%
    inner_join(p53[,c('eid','followup_time')], by = "eid") %>%
    inner_join(age_v0, by = "eid") %>%
    inner_join(protein, by = "eid")
  names(data)
  
  # follow-up time 转成年
  data <- data %>% mutate(followup_year = followup_time / 365.25)
  
  # 校正变量
  covariates <- c("age_i0", "followup_year", "sex", "Mixed", "Asian", "Black", "Chinese", "others", "education_class")
  
  protein_cols <- setdiff(names(protein), "eid")
  
  results <- lapply(seq_along(protein_cols), function(i) {
    
    protein_name <- protein_cols[i]
    
    # 当前分析所需变量
    vars <- c("age_gap_adjusted", protein_name, covariates)
    
    d <- data %>%
      dplyr::select(all_of(vars)) %>%
      filter(complete.cases(.))
    
    # 回归模型
    formula <- as.formula(paste("age_gap_adjusted ~", paste(c(protein_name, covariates), collapse = " + ")))
    fit <- lm(formula, data = d)
    
    # 提取 protein 的结果
    coef_table <- summary(fit)$coefficients
    x <- coef_table[protein_name, ]
    
    data.frame(
      protein = protein_name,
      N = nrow(d),
      beta = x["Estimate"],
      SE = x["Std. Error"],
      t = x["t value"],
      P = x["Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
    
  })
  
  results <- bind_rows(results)
  results <- results %>% mutate(FDR = p.adjust(P, method = "BH"))
  
  # 7. 保存结果
  outcome_name <- tools::file_path_sans_ext(basename(csv_file))
  outcome_name <- str_remove(outcome_name,"_oof$")
  
  results$outcome <- outcome_name
  all_results[[outcome_name]] <- results
  
  fwrite(results, paste0("result/F6_蛋白代谢相关/protein_", outcome_name,".csv"))
}






# 代谢====

# 所有代谢数据
meta = fread('/path/to/shared_project/data/res_data/all_meta_scaled.csv')
names(meta)

# 随访时间 V0-V2
p53 <- fread('shared_project/data/p53.csv')
p53 <- p53 %>%
  mutate(across(where(is.character), ~ if_else(str_detect(., "^\\s*$"), NA_character_, .))) %>%
  mutate(followup_time = as.numeric(as.Date(p53_i2) - as.Date(p53_i0)))

# V0 age
age_v0 <- fread('shared_project/data/sex_age_race.csv')
age_v0 <- age_v0[,c(1,3)]
# 预测年龄
csv_files <- list.files(path = 'data', 
                        pattern = ".*oof.*\\.csv$", full.names = TRUE)

all_results <- list()

for (csv_file in csv_files) {
  pred_age <- fread(csv_file)
  
  data <- pred_age %>%
    inner_join(p53[,c('eid','followup_time')], by = "eid") %>%
    inner_join(age_v0, by = "eid") %>%
    inner_join(meta, by = "eid")
  names(data)
  
  # follow-up time 转成年
  data <- data %>% mutate(followup_year = followup_time / 365.25)
  
  # 校正变量
  covariates <- c("age_i0", "followup_year", "sex", "Mixed", "Asian", "Black", "Chinese", "others", "education_class")
  
  meta_cols <- setdiff(names(meta), "eid")
  
  results <- lapply(seq_along(meta_cols), function(i) {
    
    meta_name <- meta_cols[i]
    
    # 当前分析所需变量
    vars <- c("age_gap_adjusted", meta_name, covariates)
    
    d <- data %>%
      dplyr::select(all_of(vars)) %>%
      filter(complete.cases(.))
    
    # 回归模型
    formula <- as.formula(paste("age_gap_adjusted ~", paste(c(meta_name, covariates), collapse = " + ")))
    print(formula)
    fit <- lm(formula, data = d)
    
    # 提取 meta 的结果
    coef_table <- summary(fit)$coefficients
    x <- coef_table[meta_name, ]
    
    data.frame(
      meta = meta_name,
      N = nrow(d),
      beta = x["Estimate"],
      SE = x["Std. Error"],
      t = x["t value"],
      P = x["Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
    
  })
  
  results <- bind_rows(results)
  results <- results %>% mutate(FDR = p.adjust(P, method = "BH"))
  
  # 7. 保存结果
  outcome_name <- tools::file_path_sans_ext(basename(csv_file))
  outcome_name <- str_remove(outcome_name,"_oof$")
  
  results$outcome <- outcome_name
  all_results[[outcome_name]] <- results
  
  fwrite(results, paste0("result/F6_蛋白代谢相关/meta_", outcome_name,".csv"))
}


# 蛋白代谢火山图====
library(ggplot2)
library(dplyr)
## 蛋白====

csv_files <- list.files(path = 'result/F6_蛋白代谢相关', 
                        pattern = ".*protein.*\\.csv$", full.names = TRUE)


for (csv_file in csv_files) {
  all_pro <- fread(csv_file)
  all_pro$significant <- ifelse(all_pro$FDR<0.05,T,F)
  all_pro$lg10FDR <- -log10(all_pro$FDR)
  head(all_pro)
  
  
  sig_up <- all_pro[all_pro$significant & all_pro$beta > 0, ]
  sig_down <- all_pro[all_pro$significant & all_pro$beta < 0, ]
  top5_up <- sig_up[order(sig_up$FDR), ][1:5, ]
  top5_down <- sig_down[order(sig_down$FDR), ][1:5, ]
  top10_custom <- rbind(top5_up, top5_down)
  
  
  library(ggrepel)
  # 绘图
  p <- ggplot(all_pro, aes(x = beta, y = lg10FDR)) +
    geom_point(aes(color = ifelse(significant, ifelse(beta > 0, "#e94234", "#269846"), "gray")), size = 1.8, alpha = 0.7) +
    labs(x = "Standardized coefficient", y = expression(bold(-log[10])~bold('(')~bold('FDR)'))) +
    theme_minimal() +
    scale_color_identity() +  # 直接使用手动指定的颜色
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
               color = "black", linewidth = 0.6) +
    theme_classic() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title = element_text(size = 18,family = "Times New Roman", face='bold', color='black'),  
      axis.text.x = element_text(size = 12,family = "Times New Roman", face='bold', color='black'),
      axis.text.y = element_text(size = 12,family = "Times New Roman", face='bold', color='black'),
    )  + 
    geom_text_repel(data = top10_custom, aes(label = protein), color="black", size=3.8, fontface="bold", family="Times",
                    arrow = arrow(ends="first", length = unit(0.01, "npc")), box.padding = 0.2,
                    point.padding = 0.3, segment.color = 'black', segment.size = 0.3, force = 1, max.iter = 3e3)
  
  ggsave(gsub('csv', 'png', csv_file), p, width = 6.2,height = 5,dpi = 900)
}
## 代谢====
csv_files <- list.files(path = 'result/F6_蛋白代谢相关', 
                        pattern = ".*meta.*\\.csv$", full.names = TRUE)


for (csv_file in csv_files) {
  all_meta <- fread(csv_file)
  all_meta$significant <- ifelse(all_meta$FDR<0.05,T,F)
  all_meta$lg10FDR <- -log10(all_meta$FDR)
  head(all_meta)
  
  
  sig_up <- all_meta[all_meta$significant & all_meta$beta > 0, ]
  sig_down <- all_meta[all_meta$significant & all_meta$beta < 0, ]
  top5_up <- sig_up[order(sig_up$FDR), ][1:5, ]
  top5_down <- sig_down[order(sig_down$FDR), ][1:5, ]
  top10_custom <- rbind(top5_up, top5_down)
  
  
  library(ggrepel)
  # 绘图
  p <- ggplot(all_meta, aes(x = beta, y = lg10FDR)) +
    geom_point(aes(color = ifelse(significant, ifelse(beta > 0, "#e94234", "#269846"), "gray")), size = 1.8, alpha = 0.7) +
    labs(x = "Standardized coefficient", y = expression(bold(-log[10])~bold('(')~bold('FDR)'))) +
    theme_minimal() +
    scale_color_identity() +  # 直接使用手动指定的颜色
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
               color = "black", linewidth = 0.6) +
    theme_classic() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title = element_text(size = 18,family = "Times New Roman", face='bold', color='black'),  
      axis.text.x = element_text(size = 12,family = "Times New Roman", face='bold', color='black'),
      axis.text.y = element_text(size = 12,family = "Times New Roman", face='bold', color='black'),
    )  + 
    geom_text_repel(data = top10_custom, aes(label = meta), color="black", size=3.8, fontface="bold", family="Times",
                    arrow = arrow(ends="first", length = unit(0.01, "npc")), box.padding = 0.2,
                    point.padding = 0.3, segment.color = 'black', segment.size = 0.3, force = 1, max.iter = 3e3)
  
  ggsave(gsub('csv', 'png', csv_file), p, width = 6.2,height = 5,dpi = 900)
}

# 蛋白KEGG-GO ====
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(reactome.db)
library(gground)
library(ggprism)
library(tidyverse)


csv_files <- list.files(path = 'result/F6_蛋白代谢相关', 
                        pattern = ".*protein.*\\.csv$", full.names = TRUE)


for (csv_file in csv_files) {
  all_pro <- fread(csv_file)
  all_pro$significant <- ifelse(all_pro$FDR<0.05,T,F)
  all_pro$lg10FDR <- -log10(all_pro$FDR)
  head(all_pro)
  data <- all_pro %>% filter(FDR < 0.05)
  interactors <- data$protein 

  
  interactors <- unique(interactors)
  
  # 1. SYMBOL -> ENTREZID
  id <- bitr(interactors, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # 已经成功匹配的 SYMBOL
  matched <- unique(id$SYMBOL)
  
  # 2. 找出 SYMBOL 没有匹配上的基因
  fail <- setdiff(interactors, matched)
  
  length(interactors)
  length(matched)
  length(fail)
  
  # 3. 对 fail 使用 ALIAS -> ENTREZID
  geneID_add <- AnnotationDbi::select(org.Hs.eg.db, keys = fail, columns = c("ENTREZID"), keytype = "ALIAS")
  
  # 去掉没有匹配上的
  geneID_add <- geneID_add %>% filter(!is.na(ENTREZID))
  
  # 4. 合并 SYMBOL 和 ALIAS 匹配结果
  id_add <- geneID_add %>%
    dplyr::rename(SYMBOL = ALIAS) %>%
    dplyr::select(SYMBOL, ENTREZID)
  
  id_final <- bind_rows(id, id_add) %>%
    distinct(SYMBOL, ENTREZID)
  
  # go
  go <- enrichGO(gene          = id_final$ENTREZID,
                 OrgDb         = "org.Hs.eg.db",
                 ont           = "ALL",#"BP", "MF", "CC"
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)
  
  GO <- go@result 
  write.csv(GO, paste0('result/F6_富集分析/GO_', basename(csv_file)), row.names = F)
  # kegg
  R.utils::setOption("clusterProfiler.download.method","auto")
  
  kegg <- enrichKEGG(
    gene = id_final$ENTREZID,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    organism = "hsa")  
  KEGG <- kegg@result
  write.csv(KEGG, paste0('result/F6_富集分析/KEGG_', basename(csv_file)), row.names = F)
  
}

# 绘图
traits <- c('double_right', 'single_right_hip', 'single_right_knee')
path_folder <- 'result/F6_富集分析/'
for (trait in traits) {
  GO <- fread(paste0(path_folder, 'GO_protein_', trait, '.csv'))
  KEGG <- fread(paste0(path_folder, 'KEGG_protein_', trait, '.csv'))
  
  use_pathway <- GO %>%
    group_by(ONTOLOGY) %>%
    arrange(p.adjust, pvalue, desc(Count), .by_group = TRUE) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    bind_rows(
      KEGG %>%
        arrange(p.adjust, pvalue, desc(Count)) %>%
        slice_head(n = 5) %>%
        mutate(ONTOLOGY = "KEGG")) %>%
    mutate(
      ONTOLOGY = case_when(
        ONTOLOGY == "BP" ~ "GO: Biological process",
        ONTOLOGY == "CC" ~ "GO: Cellular component",
        ONTOLOGY == "MF" ~ "GO: Molecular function",
        TRUE ~ ONTOLOGY)) %>%
    mutate(
      ONTOLOGY = factor(
        ONTOLOGY,
        levels = c(
          "KEGG",
          "GO: Molecular function",
          "GO: Cellular component",
          "GO: Biological process"
        )
      ),
      neg_log10_p_adjust = -log10(p.adjust)
    ) %>%
    arrange(ONTOLOGY, neg_log10_p_adjust) %>%
    mutate(Description = factor(Description, levels = Description)) %>%
    tibble::rowid_to_column("index")
  
  # 左侧分类标签和基因数量点图的宽度
  width <- 0.3
  # x 轴长度
  xaxis_max <- max(-log10(use_pathway$p.adjust)) + 1
  # 左侧分类标签数据
  rect.data <- group_by(use_pathway, ONTOLOGY) %>%
    reframe(n = n()) %>%
    ungroup() %>%
    mutate(
      xmin = -3 * width,
      xmax = -2 * width,
      ymax = cumsum(n),
      ymin = lag(ymax, default = 0) + 0.6,
      ymax = ymax + 0.4
    )
  
  
  # 修正颜色向量：与柱状图颜色顺序一致
  pal <- c(
    'KEGG' = '#7bc4e2',            # KEGG = 蓝色
    'GO: Molecular function' = '#acd372',  # MF = 绿色
    'GO: Cellular component' = '#fbb05b',  # CC = 橙色
    'GO: Biological process' = '#F1766D'   # BP = 红色
  )
  
  count_min <- min(use_pathway$Count, na.rm = TRUE)
  count_max <- max(use_pathway$Count, na.rm = TRUE)
  
  count_breaks <- pretty(c(count_min, count_max), n = 2)
  
  count_breaks <- count_breaks[count_breaks >= count_min & count_breaks <= count_max]
  
  # 如果 pretty() 最终不足3个，就直接用最小、中间、最大
  if (length(count_breaks) < 3) {
    count_breaks <- unique(round(seq(count_min, count_max, length.out = 3)))
  }
  
  p <- ggplot(use_pathway, aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
    geom_round_col(aes(y = Description), width = 0.6, alpha = 0.8) +
    geom_text(aes(x = 0.05, label = str_replace(Description, "^([a-zA-Z])", toupper)),
              hjust = 0, size = 9, fontface = "bold", family = "Times New Roman") +
    geom_point(aes(x = -width, size = Count), shape = 21) +
    geom_text(aes(x = -width, label = Count), size = 7, 
              fontface = "bold", family = "Times New Roman") +
    scale_size_continuous(
      name = "Gene count", range = c(10, 20), 
      breaks = count_breaks, limits = c(count_min, count_max)) +
    geom_text(aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = ''),
              data = rect.data,fontface = "bold",family = "Times New Roman",
              inherit.aes = FALSE) +
    geom_segment(aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
                 linewidth = 1.5, inherit.aes = FALSE) +
    labs(y = NULL, x = expression(bold(-log[10])~bold('(')~bolditalic('P')~bold('adjust)'))) +
    scale_fill_manual(name = 'Category', values = pal) +
    scale_x_continuous(breaks = seq(0, xaxis_max, 2), expand = expansion(c(0, 0))) +
    theme_prism() +
    theme(
      text = element_text(family = "Times New Roman", face = 'bold'),
      axis.text = element_text(family = "Times New Roman", size = 22),
      axis.title = element_text(family = "Times New Roman", size = 26),
      legend.title = element_text(family = "Times New Roman", size = 24),
      legend.text = element_text(family = "Times New Roman", face = 'bold', size = 22),
      axis.line = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()) +
    guides(fill = guide_legend(reverse = TRUE, keyheight = 2,  keywidth = 2, 
                               label.position = "right", label.vjust = 0.5))
  
  ggsave(paste0(path_folder, 'protein_', trait, '.pdf'),p, width=16,
         height=14, device = cairo_pdf)
}



