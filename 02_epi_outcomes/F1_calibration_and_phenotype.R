# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 02_epi_outcomes
# 核心：预测年龄线性校准 → 计算年龄加速 age_gap → 合并结局与 BMD，生成下游主键表
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
library(data.table)
library(dplyr)
library(survival)
library(ggplot2)

# 校正年龄 ====

# 疾病数据
disease_data <- fread("data/incident_death_disease_v2.csv")

# 协变量
cov_data <- fread("data/edu_race.csv")

# 骨密度
bmd_data <- fread('data/DXA_data_i2.csv')
bmd_data <- bmd_data[,c('eid', 'L1_L4_BMD', 'Femur_neck_BMD_right')]
eid <- fread('/path/to/shared_project/data/eid_241031.csv')
bmd_data <- inner_join(eid,bmd_data)
bmd_data <- bmd_data[,-1]
colnames(bmd_data)[1] <- 'eid'

# 预测年龄
csv_files <- list.files(path = 'data/age_prediction', 
                        pattern = "\\.csv$", full.names = TRUE)
csv_files
# 循环cox
result_list <- list()
for (csv_file in csv_files) {
  pred_age <- fread(csv_file)
  cat("Prediction age rows:", nrow(pred_age), "\n")
  
  pred_age$eid <- as.character(pred_age$eid)
  disease_data$eid <- as.character(disease_data$eid)
  cov_data$eid <- as.character(cov_data$eid)
  bmd_data$eid <- as.character(bmd_data$eid)
  
  # 校验是否缺失变量
  required_pred_age <- c("eid", "age", "predicted_age")
  missing_pred_age <- setdiff(required_pred_age, names(pred_age))
  
  if (length(missing_pred_age) > 0) {
    stop("\nMissing columns in pred_age:\n",
         paste(missing_pred_age, collapse = ", "))
  }
  
  if (anyDuplicated(pred_age$eid) > 0) {
    stop("\nDuplicate eid detected in pred_age!\n")
  }
  
  
  # 预测年龄校正
  fit <- lm(predicted_age ~ age, data = pred_age)
  
  alpha <- coef(fit)[["age"]]
  beta  <- coef(fit)[["(Intercept)"]]
  
  pred_age <- pred_age %>%
    mutate(
      adjusted_age = predicted_age + (age - (alpha * age + beta)),
      age_gap_adjusted = adjusted_age - age,
      age_i2 = age
    )
  
  
  data <- pred_age %>%
    inner_join(disease_data, by = "eid") %>%
    inner_join(cov_data, by = "eid") %>%
    inner_join(bmd_data, by = "eid")
  
  cat("\nFinal merged N:", nrow(data), "\n")
  
  
  write.csv(data, paste0('data/', gsub('data/age_prediction/', '', csv_file)),row.names = F)
}


# 预测年龄相关性点图====
# 预测年龄
library(ggdensity)
library(ggblanket)
library(ggsci)
csv_files <- list.files(path = 'data', 
                        pattern = ".*oof.*\\.csv$", full.names = TRUE)

for (csv_file in csv_files) {
  
  pred_age <- fread(csv_file)
  
  pred_age <- pred_age %>%
    mutate(
      sex_group = case_when(
        sex == 0 ~ "Female",
        sex == 1 ~ "Male",
        TRUE ~ as.character(sex)
      ),
      model = basename(csv_file)
    ) %>%
    dplyr::select(eid, age, predicted_age, adjusted_age, sex_group) 
  
  # 添加 Total
  pred_age_plot <- bind_rows(pred_age %>% mutate(sex_group = "Total"), pred_age) %>%
    mutate(group = factor(sex_group, levels = c("Total", "Female", "Male")))
  
  # 计算各组样本量
  n_total  <- sum(pred_age_plot$sex_group == "Total")
  n_female <- sum(pred_age_plot$sex_group == "Female")
  n_male   <- sum(pred_age_plot$sex_group == "Male")
  
  # 修改 facet 标签
  facet_labels <- c(
    "Total"  = paste0("Total (n=", n_total, ")"),
    "Female" = paste0("Female (n=", n_female, ")"),
    "Male"   = paste0("Male (n=", n_male, ")")
  )
  
  
  p <- ggplot(pred_age_plot, aes(x = age, y = adjusted_age, fill = group)) +
    geom_hdr(probs = c(0.99, 0.95,0.9, 0.8,0.7,0.6, 0.5,0.4,0.3,0.2,0.1)) + 
    scale_fill_aaas() +
    facet_wrap(vars(group), labeller = as_labeller(facet_labels)) +
    theme_bw() +
    scale_x_continuous(limits = c(40, 90)) +
    scale_y_continuous(limits = c(40, 90)) +
    labs(title=NULL, x="Chronological age (years)", y="Predicted age (years)") +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title = element_text(size = 18, family = "Times New Roman", face = 'bold', color = 'black'),
      axis.text.y = element_text(size = 14, family = "Times New Roman", face = 'bold', color = 'black'),
      axis.text.x = element_text(size = 14, family = "Times New Roman", face = 'bold', color = 'black'),
      panel.border = element_rect(size = 1.5, color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1, "lines"),
      strip.text = element_text(size = 16, family = "Times New Roman", face = 'bold', color = 'black'),
      strip.background = element_rect(color = "black", size = 1.5),
      legend.position = 'none'
    )
  
  ggsave(paste0('result/F1_密度图/',gsub('csv', 'png', basename(csv_file))), p, width = 12,height = 5,dpi = 600)
}


# 年龄加速相关性 ====

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

library(Hmisc)
library(dplyr)

double_right <- fread("data/double_right_oof.csv") %>%
  dplyr::select(eid, age_gap_adjusted) %>%
  dplyr::rename(double_right = age_gap_adjusted)

right_knee <- fread("data/single_right_knee_oof.csv") %>%
  dplyr::select(eid, age_gap_adjusted) %>%
  dplyr::rename(single_right_knee = age_gap_adjusted)

right_hip <- fread("data/single_right_hip_oof.csv") %>%
  dplyr::select(eid, age_gap_adjusted) %>%
  dplyr::rename(single_right_hip = age_gap_adjusted)

double_left <- fread("data/double_left_external_valid.csv") %>%
  dplyr::select(eid, age_gap_adjusted) %>%
  dplyr::rename(double_left = age_gap_adjusted)

left_knee <- fread("data/single_left_knee_external_valid.csv") %>%
  dplyr::select(eid, age_gap_adjusted) %>%
  dplyr::rename(single_left_knee = age_gap_adjusted)

left_hip <- fread("data/single_left_hip_external_valid.csv") %>%
  dplyr::select(eid, age_gap_adjusted) %>%
  dplyr::rename(single_left_hip = age_gap_adjusted)


names(double_left)




age_acc <- double_right %>%
  inner_join(right_knee, by = "eid") %>%
  inner_join(right_hip, by = "eid") %>%
  left_join(double_left, by = "eid") %>%
  left_join(left_knee, by = "eid") %>%
  left_join(left_hip, by = "eid")
colnames(age_acc)

vars <- c("double_right", "single_right_knee", "single_right_hip",
  "double_left", "single_left_knee", "single_left_hip")


# 1. Correlation
cor_result <- rcorr(as.matrix(
  age_acc %>% dplyr::select(all_of(vars))), type = "pearson")

# 3. 两两比较
pairs <- combn(vars, 2, simplify = FALSE)

result_table <- lapply(pairs, function(x) {
  v1 <- x[1]
  v2 <- x[2]

    # 只保留两个变量都有数据的样本
  d <- age_acc %>%
    dplyr::select(all_of(c(v1, v2))) %>%
    filter(complete.cases(.))
  
  # Mean difference
  diff <- abs(d[[v1]] - d[[v2]])
  
  data.frame(
    Comparison = paste(v1, "vs", v2),
    N = nrow(d),
    Pearson_r = cor_result$r[v1, v2],
    P_value = cor_result$P[v1, v2],
    Mean_difference = mean(diff)
  )
}) %>%
  bind_rows()


# 3. 添加自身 vs 自身
self_result <- lapply(vars, function(v) {
  
  # 该变量自身非缺失的样本
  d <- age_acc %>%
    dplyr::select(all_of(v)) %>%
    filter(!is.na(.data[[v]]))
  
  data.frame(
    Comparison = paste(v, "vs", v),
    N = nrow(d),
    Pearson_r = 1,
    P_value = NA_real_,
    Mean_difference = 0
  )
}) %>%
  bind_rows()

# 4. 合并
result_table <- bind_rows(result_table, self_result) %>%
  mutate(
    Pearson_r = round(Pearson_r, 8),
    Mean_difference = round(Mean_difference, 4),
    P = P_value,
    P_value = ifelse(
      is.na(P_value),
      NA_character_,
      format.pval(P_value, digits = 3, eps = 0.001)
    )
  )

result_table



write.csv(result_table, 'result/F1_密度图/cor.csv',row.names = F)

result_table <- fread('result/F1_密度图/cor.csv')
result_table


library(data.table)
library(dplyr)
library(ggplot2)
library(grid)

# 1. 整理 result_table
result <- result_table %>%
  mutate(
    var1 = sub(" vs .*", "", Comparison),
    var2 = sub(".* vs ", "", Comparison),
    r = Pearson_r,
    MAD = Mean_difference
  )

result

# A. 三种年龄加速相互间相关性
A_vars <- c("single_right_knee", "single_right_hip", "double_right")

A_data <- result %>%
  filter(var1 %in% A_vars, var2 %in% A_vars) %>%
  mutate(
    Var1 = factor(var1, levels = A_vars),
    Var2 = factor(var2, levels = A_vars),
    label = sprintf("%.3f", MAD)
  )

my_colors <- c("#fbe6e6", "#f8c8c8", "#d45959" )
pA <- ggplot(A_data, aes(x = Var2, y = Var1, fill = r)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 6, lineheight = 0.9, 
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
        size = 18, face = "bold", colour = "black"
      ),
      barwidth  = unit(0.8, "cm"),
      barheight = unit(10, "cm")      # ← 色条加长
    )
  ) +
  # X轴：从左到右
  scale_x_discrete(limits = A_vars, labels = c('Knee', 'Hip', 'Knee - Hip')) +
  # Y轴：从上到下
  scale_y_discrete(limits = rev(A_vars),labels = c('Knee - Hip', 'Hip', 'Knee')) +
  coord_fixed() +
  theme_minimal() +
  theme(
    plot.title = element_text(family = "Times New Roman", size = 20, 
                              face = "bold", colour = "black", hjust = 0.5),
    legend.title = element_text(size = 18, family = "Times New Roman", 
                                colour = "black", face = "bold"),
    legend.text = element_text(family = "Times New Roman", size = 16, face = "bold"),
    axis.text.y = element_text(size = 18, family = "Times New Roman", colour = "black", face = "bold"),
    axis.text.x = element_text(size = 18, family = "Times New Roman", colour = "black", face = "bold"),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white", color = "white"),
    legend.key.size = unit(0.8, "cm")
  )

pA

ggsave('result/F1_密度图/right.pdf', pA, height = 7, width = 10, device = cairo_pdf)


# A. 三种年龄加速相互间相关性
B_vars <- c("single_right_knee vs single_left_knee", 
            "single_right_hip vs single_left_hip",
            "double_right vs double_left")

B_vars <- c("single_right_knee vs single_left_knee", 
            "single_right_hip vs single_left_hip",
            "double_right vs double_left")

B_data <- result %>%
  filter(Comparison %in% B_vars) %>%
  mutate(
    Var1 = factor(var1, levels = A_vars),
    Var2 = factor(var2, levels = gsub('right', 'left', A_vars)),
    label = sprintf("%.3f", MAD)) 

my_colors <- c("#d0e8ff", "#c0e0f8", "#5a97d0" )
pB <- ggplot(B_data, aes(x = Var2, y = Var1, fill = r)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 6, lineheight = 0.9, 
            family = "Times New Roman", colour = "black", fontface = "bold") +
  scale_fill_gradientn(
    colors = my_colors,
    limits = c(0, 0.6),
    breaks = c(0, 0.2, 0.4, 0.6),
    na.value = "white",
    name = "Correlation coefficients (r)",
    guide = guide_colorbar(
      title.position = "left",        # 标题放到色条左侧
      title.theme = element_text(
        angle = 90,                   # 竖向
        hjust = 0.5,                  # 在色条长度方向上居中
        vjust = 0.5,
        family = "Times New Roman",
        size = 18, face = "bold", colour = "black"
      ),
      barwidth  = unit(0.8, "cm"),
      barheight = unit(10, "cm")      # ← 色条加长
    )
  ) +
  # X轴：从左到右
  scale_x_discrete(
    limits = rev(gsub('right','left',A_vars)),
    labels = c('Knee - Hip', 'Hip', 'Knee')
  ) +
 
  # Y轴：从上到下
  scale_y_discrete(
    limits = rev(A_vars),
    labels = c('Knee - Hip', 'Hip', 'Knee')
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    plot.title = element_text(family = "Times New Roman", size = 20, 
                              face = "bold", colour = "black", hjust = 0.5),
    legend.title = element_text(size = 18, family = "Times New Roman", 
                                colour = "black", face = "bold"),
    legend.text = element_text(family = "Times New Roman", size = 16, face = "bold"),
    axis.text.y = element_text(size = 18, family = "Times New Roman", colour = "black", face = "bold"),
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1, vjust = 1, 
                               family = "Times New Roman", colour = "black", face = "bold"),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white", color = "white"),
    legend.key.size = unit(0.8, "cm")
  )

pB

ggsave('result/F1_密度图/right_left.pdf', pB, height = 7, width = 10, device = cairo_pdf)

p <- pA + pB
p

ggsave('result/F1_密度图/all_right_left.pdf', p, height = 7, width = 20, device = cairo_pdf)


C_data <- result %>%
  filter(!( (var1 %in% A_vars & var2 %in% A_vars) | (Comparison %in% B_vars) )) %>%
  mutate(
    Var1 = factor(
      var1, levels = c('double_right', 'single_right_knee', 
                       'single_right_hip',  'double_left', 'single_left_knee')),
    Var2 = factor(var2, levels = c('double_left', 
                                   'single_left_knee', 'single_left_hip')),   # 或使用 unique 值
    label = sprintf("%.3f", MAD))

my_colors <- c("#d7e9bc", "#aac576", "#507d39" )
pC <- ggplot(C_data, aes(x = Var2, y = Var1, fill = r)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 6, lineheight = 0.9, 
            family = "Times New Roman", colour = "black", fontface = "bold") +
  scale_fill_gradientn(
    colors = my_colors,
    limits = c(0, 0.8),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8),
    na.value = "white",
    name = "Correlation coefficients (r)"
  ) +
  # X轴：从左到右
  scale_x_discrete(
    limits = c('double_left', 'single_left_knee', 'single_left_hip'),
    labels = c("Double left", "Single left knee", "Single left hip")) +
  # Y轴：从上到下
  scale_y_discrete(
    limits = c('double_right', 'single_right_knee', 
               'single_right_hip',  'double_left', 'single_left_knee'),
    labels = c("Double right", "Single right knee", "Single right hip",
               "Double left", "Single left knee")) +
  coord_fixed() +
  theme_minimal() +
  theme(
    plot.title = element_text(family = "Times New Roman", size = 20, 
                              face = "bold", colour = "black", hjust = 0.5),
    legend.title = element_text(size = 18, family = "Times New Roman", 
                                colour = "black", face = "bold"),
    legend.text = element_text(family = "Times New Roman", size = 16, face = "bold"),
    axis.text.y = element_text(size = 18, family = "Times New Roman", colour = "black", face = "bold"),
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1, vjust = 1, 
                               family = "Times New Roman", colour = "black", face = "bold"),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white", color = "white"),
    legend.key.size = unit(0.8, "cm")
  )

pC

ggsave('result/F1_密度图/right_left_other.pdf', pC, height = 7, width = 10, device = cairo_pdf)


# cox 未校正骨密度 右侧====

# 疾病映射
disease_mapping <- c(
  "death" = "Death",
  "CVD" = "Cardiovascular diseases",
  "respiratory" = "Respiratory diseases",
  "digestive" = "Digestive diseases",
  "endocrine" = "Endocrine diseases",
  "neurological" = "Neurological diseases",
  "psychiatric" = "Psychiatric diseases",
  "musculoskeleta" = "Musculoskeletal diseases"
)

# 颜色
color_palette <- c("#DC143C","#1E90FF","#05004e","#007b43","#f8b500",
                   "#8B008B","#dbd0e6","#FF69B4","#ff847c","#a0a003",
                   "#08d9d6","#ba52ed","#323232")
#  Cox协变量
base_covariates <- c("age_i2", "sex", "Mixed", "Asian", "Black", 
                     "Chinese", "others", "education_class")

# 预测年龄
csv_files <- list.files(path = 'data', 
                        pattern = ".*oof.*\\.csv$", full.names = TRUE)
csv_files
# 循环cox
result_list <- list()
for (csv_file in csv_files) {
  data <- fread(csv_file)
  cat("Prediction age rows:", nrow(data), "\n")
  
  # 疾病列表
  disease_names <- c(
    "death", "CVD", "respiratory", "digestive", "endocrine", "neurological", 
    "psychiatric", "musculoskeleta"
  )
  for (disease in disease_names) {
    
    exposure <- 'age_gap_adjusted'
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
      mode = basename(csv_file),
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

write.csv(result_list, 'result/F1_cox/all.cox.csv', row.names = F)

result_list <- fread('result/F1_cox/all.cox.csv')
plot_data <- result_list %>%
  filter(grepl('right',mode)) %>%
  mutate(
    mode = case_when(
      grepl("^double_right", mode) ~ "Right hip and knee",
      grepl("^single_right_hip", mode) ~ "Right hip",
      grepl("^single_right_knee", mode) ~ "Right knee",
      TRUE ~ mode
    ),
    outcome = factor(
      outcome, 
      levels = c("death", "CVD", "respiratory", "digestive", "endocrine", 
                 "neurological", "psychiatric", "musculoskeleta")
    )
  ) %>%
  mutate(
    mode = factor(mode, levels = c("Right knee", "Right hip", "Right hip and knee")),
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
    mode = dplyr::recode(
      mode,
      `Right knee` = "Knee",
      `Right hip` = "Hip",
      `Right hip and knee` = "Knee - Hip"
    )
  )
p <- ggplot(plot_data, aes(x = outcome, y = HR, color = mode)) +
  geom_hline(yintercept = 1,colour = "grey20",linetype = "dashed",linewidth = 0.6) +
  geom_linerange(aes(ymin = CI_low, ymax = CI_high),
                 position = position_dodge(0.6),linewidth = 2) +
  geom_point(position = position_dodge(0.6), shape = 15, size = 5 ) +
  scale_y_continuous(limits = c(0.99, 1.08), breaks = seq(1, 1.08, 0.02)) +
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

ggsave('result/F1_cox/right.png', p, width = 9, height = 7,  dpi = 600)


# cox 未校正骨密度 左侧====

# 疾病映射
disease_mapping <- c(
  "death" = "Death",
  "CVD" = "Cardiovascular diseases",
  "respiratory" = "Respiratory diseases",
  "digestive" = "Digestive diseases",
  "endocrine" = "Endocrine diseases",
  "neurological" = "Neurological diseases",
  "psychiatric" = "Psychiatric diseases",
  "musculoskeleta" = "Musculoskeletal diseases"
)

# 颜色
color_palette <- c("#DC143C","#1E90FF","#05004e","#007b43","#f8b500",
                   "#8B008B","#dbd0e6","#FF69B4","#ff847c","#a0a003",
                   "#08d9d6","#ba52ed","#323232")
#  Cox协变量
base_covariates <- c("age_i2", "sex", "Mixed", "Asian", "Black", 
                     "Chinese", "others", "education_class")

# 预测年龄
csv_files <- list.files(path = 'data', 
                        pattern = ".*_external_valid.*\\.csv$", full.names = TRUE)
csv_files
# 循环cox
result_list <- list()
for (csv_file in csv_files) {
  data <- fread(csv_file)
  cat("Prediction age rows:", nrow(data), "\n")
  
  # 疾病列表
  disease_names <- c(
    "death", "CVD", "respiratory", "digestive", "endocrine", "neurological", 
    "psychiatric", "musculoskeleta"
  )
  for (disease in disease_names) {
    
    exposure <- 'age_gap_adjusted'
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
      mode = basename(csv_file),
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

write.csv(result_list, 'result/F1_cox/all_left.cox.csv', row.names = F)

result_list <- fread('result/F1_cox/all_left.cox.csv')
plot_data <- result_list %>%
  filter(grepl('left',mode)) %>%
  mutate(
    mode = case_when(
      grepl("^double_left", mode) ~ "Left hip and knee",
      grepl("^single_left_hip", mode) ~ "Left hip",
      grepl("^single_left_knee", mode) ~ "Left knee",
      TRUE ~ mode
    ),
    outcome = factor(
      outcome, 
      levels = c("death", "CVD", "respiratory", "digestive", "endocrine", 
                 "neurological", "psychiatric", "musculoskeleta")
    )
  ) %>%
  mutate(
    mode = factor(mode, levels = c("Left knee", "Left hip", "Left hip and knee")),
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
  )  %>%
  mutate(
    mode = dplyr::recode(
      mode,
      `Left knee` = "Knee",
      `Left hip` = "Hip",
      `Left hip and knee` = "Knee - Hip"
    )
  )

p <- ggplot(plot_data, aes(x = outcome, y = HR, color = mode)) +
  geom_hline(yintercept = 1,colour = "grey20",linetype = "dashed",linewidth = 0.6) +
  geom_linerange(aes(ymin = CI_low, ymax = CI_high),
                 position = position_dodge(0.6),linewidth = 2.5) +
  geom_point(position = position_dodge(0.6), shape = 15, size = 5.5 ) +
  scale_y_continuous(limits = c(0.99, 1.08), breaks = seq(1, 1.08, 0.02)) +
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

ggsave('result/F1_cox/left.png', p, width = 9, height = 7,  dpi = 600)


# cox 校正Femur_neck_BMD_right====

# 疾病映射
disease_mapping <- c(
  "death" = "Death",
  "CVD" = "Cardiovascular diseases",
  "respiratory" = "Respiratory diseases",
  "digestive" = "Digestive diseases",
  "endocrine" = "Endocrine diseases",
  "neurological" = "Neurological diseases",
  "psychiatric" = "Psychiatric diseases",
  "musculoskeleta" = "Musculoskeletal diseases"
)

# 颜色
color_palette <- c("#DC143C","#1E90FF","#05004e","#007b43","#f8b500",
                   "#8B008B","#dbd0e6","#FF69B4","#ff847c","#a0a003",
                   "#08d9d6","#ba52ed","#323232")
#  Cox协变量
base_covariates <- c("age_i2", "sex", "Mixed", "Asian", "Black", 
                     "Chinese", "others", "education_class", "Femur_neck_BMD_right")

csv_files <- list.files(path = 'data', 
                        pattern = ".*oof.*\\.csv$", full.names = TRUE)
csv_files
# 循环cox
result_list <- list()
for (csv_file in csv_files) {
  data <- fread(csv_file)
  cat("Prediction age rows:", nrow(data), "\n")
  
  # 疾病列表
  disease_names <- c(
    "death", "CVD", "respiratory", "digestive", "endocrine", "neurological", 
    "psychiatric", "musculoskeleta"
  )
  for (disease in disease_names) {
    
    exposure <- 'age_gap_adjusted'
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
      mode = basename(csv_file),
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

write.csv(result_list, 'result/F1_cox/all_Femur_neck_BMD_right.cox.csv', row.names = F)

result_list <- fread('result/F1_cox/all_Femur_neck_BMD_right.cox.csv')
plot_data <- result_list %>%
  filter(grepl('right',mode)) %>%
  mutate(
    mode = case_when(
      grepl("^double_right", mode) ~ "Right hip and knee",
      grepl("^single_right_hip", mode) ~ "Right hip",
      grepl("^single_right_knee", mode) ~ "Right knee",
      TRUE ~ mode
    ),
    outcome = factor(
      outcome, 
      levels = c("death", "CVD", "respiratory", "digestive", "endocrine", 
                 "neurological", "psychiatric", "musculoskeleta")
    )
  ) %>%
  mutate(
    mode = factor(mode, levels = c("Right knee", "Right hip", "Right hip and knee")),
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
    mode = dplyr::recode(
      mode,
      `Right knee` = "Knee",
      `Right hip` = "Hip",
      `Right hip and knee` = "Knee - Hip"
    )
  )

p <- ggplot(plot_data, aes(x = outcome, y = HR, color = mode)) +
  geom_hline(yintercept = 1,colour = "grey20",linetype = "dashed",linewidth = 0.6) +
  geom_linerange(aes(ymin = CI_low, ymax = CI_high),
                 position = position_dodge(0.6),linewidth = 2.5) +
  geom_point(position = position_dodge(0.6), shape = 15, size = 5.5 ) +
  scale_y_continuous(limits = c(0.99, 1.08), breaks = seq(1, 1.08, 0.02)) +
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

ggsave('result/F1_cox/all_Femur_neck_BMD.png', p, width = 9, height = 7,  dpi = 600)


# cox 校正L1_L4_BMD====

# 疾病映射
disease_mapping <- c(
  "death" = "Death",
  "CVD" = "Cardiovascular diseases",
  "respiratory" = "Respiratory diseases",
  "digestive" = "Digestive diseases",
  "endocrine" = "Endocrine diseases",
  "neurological" = "Neurological diseases",
  "psychiatric" = "Psychiatric diseases",
  "musculoskeleta" = "Musculoskeletal diseases"
)

# 颜色
color_palette <- c("#DC143C","#1E90FF","#05004e","#007b43","#f8b500",
                   "#8B008B","#dbd0e6","#FF69B4","#ff847c","#a0a003",
                   "#08d9d6","#ba52ed","#323232")
#  Cox协变量
base_covariates <- c("age_i2", "sex", "Mixed", "Asian", "Black", 
                     "Chinese", "others", "education_class", "L1_L4_BMD")

csv_files <- list.files(path = 'data', 
                        pattern = ".*oof.*\\.csv$", full.names = TRUE)
csv_files
# 循环cox
result_list <- list()
for (csv_file in csv_files) {
  data <- fread(csv_file)
  cat("Prediction age rows:", nrow(data), "\n")
  
  # 疾病列表
  disease_names <- c(
    "death", "CVD", "respiratory", "digestive", "endocrine", "neurological", 
    "psychiatric", "musculoskeleta"
  )
  for (disease in disease_names) {
    
    exposure <- 'age_gap_adjusted'
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
      mode = basename(csv_file),
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

write.csv(result_list, 'result/F1_cox/all_L1_L4_BMD.cox.csv', row.names = F)

result_list <- fread('result/F1_cox/all_L1_L4_BMD.cox.csv')
plot_data <- result_list %>%
  filter(grepl('right',mode)) %>%
  mutate(
    mode = case_when(
      grepl("^double_right", mode) ~ "Right hip and knee",
      grepl("^single_right_hip", mode) ~ "Right hip",
      grepl("^single_right_knee", mode) ~ "Right knee",
      TRUE ~ mode
    ),
    outcome = factor(
      outcome, 
      levels = c("death", "CVD", "respiratory", "digestive", "endocrine", 
                 "neurological", "psychiatric", "musculoskeleta")
    )
  ) %>%
  mutate(
    mode = factor(mode, levels = c("Right knee", "Right hip", "Right hip and knee")),
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
    mode = dplyr::recode(
      mode,
      `Right knee` = "Knee",
      `Right hip` = "Hip",
      `Right hip and knee` = "Knee - Hip"
    )
  )
p <- ggplot(plot_data, aes(x = outcome, y = HR, color = mode)) +
  geom_hline(yintercept = 1,colour = "grey20",linetype = "dashed",linewidth = 0.6) +
  geom_linerange(aes(ymin = CI_low, ymax = CI_high),
                 position = position_dodge(0.6),linewidth = 2.5) +
  geom_point(position = position_dodge(0.6), shape = 15, size = 5.5 ) +
  scale_y_continuous(limits = c(0.99, 1.08), breaks = seq(1, 1.08, 0.02)) +
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

ggsave('result/F1_cox/all_L1_L4_BMD.png', p, width = 9, height = 7,  dpi = 600)


# RCS 右侧 ====
library(rms)
library(survival)
library(dplyr)
library(ggplot2)
library(data.table)
library(patchwork)

csv_files <- list.files(path = 'data', 
                        pattern = ".*oof.*\\.csv$", full.names = TRUE)


for (csv_file in csv_files) {
  disease_names <- c("death", "CVD", "respiratory", "digestive", "endocrine", 
                     "neurological", "psychiatric", "musculoskeleta")
  
  disease_labels <- c(
    death = "All-cause mortality",
    CVD = "CVD",
    respiratory = "Respiratory",
    digestive = "Digestive",
    endocrine = "Endocrine",
    neurological = "Neurological",
    psychiatric = "Psychiatric",
    musculoskeleta = "Musculoskeletal"
  )
  
  
  # knots 设置
  
  knot_probs <- list(
    `3 knots` = c(0.05, 0.50, 0.95),
    `4 knots` = c(0.05, 0.35, 0.65, 0.95),
    `5 knots` = c(0.05, 0.25, 0.50, 0.75, 0.95),
    `6 knots` = c(0.05, 0.20, 0.40, 0.60, 0.80, 0.95)
  )
  
  
  # 计算 knots
  
  pred_age <- fread(csv_file)
  
  x <- pred_age$age_gap_adjusted
  
  knots_list <- lapply(knot_probs, function(p) {
    quantile(x, probs = p, na.rm = TRUE, names = FALSE)
  })
  
  knots_list
  
  # reference 和绘图范围
  
  reference_value <- median(pred_age$age_gap_adjusted, na.rm = TRUE)
  
  x_range <- quantile(pred_age$age_gap_adjusted,
                      probs = c(0.01, 0.99),
                      na.rm = TRUE)
  
  x_min <- x_range[1]
  x_max <- x_range[2]
  
  cat("Reference =", reference_value, "\n")
  cat("X-axis =", x_min, "to", x_max, "\n")
  
  
  # 输出文件夹
  
  output_dir <- "result/F1_RCS"
  
  # 设置 rms datadist
  dd <- datadist(pred_age)
  dd$limits["Adjust to", "age_gap_adjusted"] <- reference_value
  options(datadist = "dd")
  
  # 保存所有结果
  all_results <- list()
  all_predictions <- list()
  
  # 循环 knots
  for (knot_name in names(knots_list)) {
    cat("\n====================================\n")
    cat(knot_name, "\n")
    cat("====================================\n")
    
    knots <- knots_list[[knot_name]]
    
    plot_list <- list()
    
    # 10. 循环疾病
    for (disease in disease_names) {
      cat("\nRunning:", disease, "\n")
      
      # 构造 Surv
      time_var <- paste0(disease, "_time")
      event_var <- disease
      
      # 分析数据
      vars_needed <- c(
        time_var, event_var, "age_gap_adjusted", "age_i2", "sex", "Mixed", "Asian", 
        "Black", "Chinese", "others", "education_class"
      )
      
      dat <- pred_age %>%
        dplyr::select(all_of(vars_needed)) %>%
        na.omit()
      
      # 构造 formula
      fml <- as.formula(
        paste0(
          "Surv(", time_var, ", ", event_var, ") ~ ",
          "rcs(age_gap_adjusted, ", length(knots), ") + ", 
          "age_i2 + sex + Mixed + Asian + Black + Chinese + others + education_class"
        )
      )
      
      
      
      # 拟合 Cox RCS
      fit <- cph(fml, data = dat, x = TRUE, y = TRUE, surv = TRUE)
      
      # ANOVA
      a <- anova(fit)
      # AIC
      model_aic <- AIC(fit)
      # 提取 P overall / P nonlinear
      a_df <- as.data.frame(a)
      row_names <- rownames(a_df)
      
      # 找到 age_gap_adjusted 对应行
      overall_row <- which(grepl("^age_gap_adjusted$", row_names))
      nonlinear_row <- which(grepl("Nonlinear", row_names, ignore.case = TRUE))
      
      # 如果找不到，设置 NA
      p_overall <- NA_real_
      p_nonlinear <- NA_real_
      
      if (length(overall_row) > 0) {
        p_col <- grep("P|Wald", colnames(a_df), ignore.case = TRUE, value = TRUE)
        if (length(p_col) > 0) {
          p_overall <- a_df[overall_row[1], p_col[length(p_col)]]
        }
      }
      
      if (length(nonlinear_row) > 0) {
        p_col <- grep("P|Wald", colnames(a_df), ignore.case = TRUE, value = TRUE)
        if (length(p_col) > 0) {
          p_nonlinear <- a_df[nonlinear_row[1], p_col[length(p_col)]]
        }
      }
      
      # Prediction
      pred <- Predict(
        fit, age_gap_adjusted = seq(x_min, x_max, length.out = 300),
        ref.zero = TRUE, fun = exp, conf.int = 0.95
      )
      pred_df <- as.data.frame(pred)
      
      # 保存 prediction
      pred_df <- pred_df %>%
        mutate(
          outcome = disease,
          outcome_label = disease_labels[disease],
          knots = knot_name,
          reference = reference_value,
          N = nrow(dat),
          events = sum(dat[[event_var]] == 1, na.rm = TRUE),
          P_overall = p_overall,
          P_nonlinear = p_nonlinear
        )
      
      all_predictions[[paste0(knot_name, "_", disease)]] <- pred_df
      
      # P value 格式
      p_overall_text <- ifelse(
        is.na(p_overall),
        "P overall = NA",
        ifelse(
          p_overall < 0.001,
          "P overall < 0.001",
          paste0("P overall = ", formatC(
            p_overall, format = "f", digits = 3
          ))
        )
      )
      
      p_nonlinear_text <- ifelse(
        is.na(p_nonlinear),
        "P nonlinear = NA",
        ifelse(
          p_nonlinear < 0.001,
          "P nonlinear < 0.001",
          paste0("P nonlinear = ", formatC(
            p_nonlinear, format = "f", digits = 3
          ))
        )
      )
      
      # 单个 disease 图
      p <- ggplot(pred_df, aes(x = age_gap_adjusted, y = yhat)) +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.20) +
        geom_line(linewidth = 0.9) +
        geom_hline(yintercept = 1,
                   linetype = "dashed",
                   linewidth = 0.5) +
        geom_vline(xintercept = reference_value,
                   linetype = "dotted",
                   linewidth = 0.5) +
        labs(
          title = disease_labels[disease],
          subtitle = paste(p_overall_text, "\n", p_nonlinear_text),
          x = "Age gap (adjusted)",
          y = "Hazard ratio (95% CI)"
        ) +
        coord_cartesian(xlim = c(x_min, x_max)) +
        theme_classic(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9),
          axis.title = element_text(face = "bold")
        )
      plot_list[[disease]] <- p
      
      # 保存统计结果
      all_results[[paste0(knot_name, "_", disease)]] <-
        data.frame(
          outcome = disease,
          outcome_label = disease_labels[disease],
          knots = knot_name,
          aic = model_aic,
          knot_number = length(knots),
          knot_positions = paste(round(knots, 4), collapse = ", "),
          N = nrow(dat),
          events = sum(dat[[event_var]] == 1, na.rm = TRUE),
          reference = reference_value,
          x_min = x_min,
          x_max = x_max,
          P_overall = p_overall,
          P_nonlinear = p_nonlinear
        )
    }
    
    # 8个疾病合并成一张图
    combined_plot <- wrap_plots(plot_list, ncol = 2) +
      plot_annotation(
        title = paste0("Restricted cubic spline analysis: ", knot_name),
        subtitle = paste0("Reference = median age gap (", round(reference_value, 3), ")")
      )
    
    # 保存图片
    file_name <- paste0(gsub("oof.csv", "", basename(csv_file)), gsub(" ", "_", knot_name))
    
    ggsave(
      filename = file.path(output_dir, paste0("RCS_", file_name, ".png")),
      plot = combined_plot, width = 7, height = 14, dpi = 300
    )
    
    # 同时保存 PDF
    ggsave(
      filename = file.path(output_dir, paste0("RCS_", file_name, ".pdf")),
      plot = combined_plot, width = 7, height = 14
    )
    cat("\nSaved:", file.path(output_dir, paste0("RCS_", file_name, ".png")), "\n")
  }
  
  # 13. 合并统计结果
  results_df <- bind_rows(all_results)
  results_df <- results_df %>%
    arrange(knot_number, match(outcome, disease_names))
  
  
  # 保存 summary
  write.csv(results_df,
            file.path(output_dir, file_name <- paste0(gsub("oof.csv", "", basename(csv_file)), "RCS_summary.csv")),
            row.names = FALSE)
  
  
  # 合并所有 prediction
  prediction_df <- bind_rows(all_predictions)
  write.csv(prediction_df,
            file.path(output_dir, file_name <- paste0(gsub("oof.csv", "", basename(csv_file)), "RCS_predictions.csv")),
            row.names = FALSE)
  
  # 查看结果
  print(results_df)
}

## 绘图 death ====
double_rcs_data <- fread('result/F1_RCS/double_right_RCS_predictions.csv')
hip_rcs_data <- fread('result/F1_RCS/single_right_hip_RCS_predictions.csv')
knee_rcs_data <- fread('result/F1_RCS/single_right_knee_RCS_predictions.csv')
head(double_rcs_data)

double_death <- double_rcs_data[outcome == "death"]
hip_death    <- hip_rcs_data[outcome == "death"]
knee_death   <- knee_rcs_data[outcome == "death"]

# 添加 trait
double_death$trait <- 'Knee - Hip'
hip_death$trait <- 'Hip'
knee_death$trait <- 'Knee'

rcs_death <- rbindlist(list(double_death, hip_death, knee_death),
                       use.names = TRUE,
                       fill = TRUE)

rcs_death$trait <- factor(rcs_death$trait, levels = c("Knee", "Hip", "Knee - Hip"))

rcs_death <- rcs_death[rcs_death$knots == '3 knots',]

# 统一坐标轴范围
y_range <- range(c(rcs_death$lower, rcs_death$upper), na.rm = TRUE)

# 如果希望上下留一点空白
y_range <- c(floor(y_range[1] * 10) / 10, ceiling(y_range[2] * 10) / 10)

# 统一主题
my_theme <- theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 26, family = "Times New Roman", face = "bold", colour = "black"),
    axis.title = element_text(size = 24, family = "Times New Roman", face = 'bold'),
    axis.text = element_text(size = 20, family = "Times New Roman", face = 'bold', colour = 'black'),
    legend.title=element_blank(),
    legend.text = element_text(size = 18, family = "Times New Roman", face = 'bold'),
    axis.line = element_line(color  = "black", linewidth = 1.2, lineend = "square")
  )


# Knee
d1 <- rcs_death[rcs_death$trait == "Knee", ]

p1 <- ggplot(d1, aes(x = age_gap_adjusted, y = yhat)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7, color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#006affff", alpha = 0.20) +
  geom_line(color = "#006affff", linewidth = 1.3) +
  scale_x_continuous(limits = c(-10,10), breaks = seq(-10,10,5)) +
  scale_y_continuous(limits = y_range) +
  labs(title = "Knee", x = "Age acceleration (years)", y = "Hazard ratio (95% CI)") +
  my_theme

p1

# Hip
d2 <- rcs_death[rcs_death$trait == "Hip", ]

p2 <- ggplot(d2, aes(x = age_gap_adjusted, y = yhat)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7, color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#1dcc2fff", alpha = 0.20) +
  geom_line(color = "#1dcc2fff", linewidth = 1.3) +
  scale_x_continuous(limits = c(-10,10), breaks = seq(-10,10,5)) +
  scale_y_continuous(limits = y_range) +
  labs(title = "Hip", x = "Age acceleration (years)", y = "Hazard ratio (95% CI)") +
  my_theme

# Knee - Hip
d3 <- rcs_death[rcs_death$trait == "Knee - Hip", ]

p3 <- ggplot(d3, aes(x = age_gap_adjusted, y = yhat)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7, color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#FD3E5Cff", alpha = 0.20) +
  geom_line(color = "#FD3E5Cff", linewidth = 1.3) +
  scale_x_continuous(limits = c(-10,10), breaks = seq(-10,10,5)) +
  scale_y_continuous(limits = y_range) +
  labs(title = "Knee - Hip", x = "Age acceleration (years)", y = "Hazard ratio (95% CI)") +
  my_theme



p <- p1 + p2 + p3 + plot_layout(ncol = 3)
p

ggsave('result/F1_RCS/rcs_right_knee.png', p1, width = 6, height = 6,  dpi = 600)
ggsave('result/F1_RCS/rcs_right_hip.png', p2, width = 6, height = 6,  dpi = 600)
ggsave('result/F1_RCS/rcs_right_double.png', p3, width = 6, height = 6,  dpi = 600)
ggsave('result/F1_RCS/rcs_right_all.png', p, width = 18, height = 6,  dpi = 600)



## 绘图 disease ====
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)

double_rcs_data <- fread('result/F1_RCS/double_right_RCS_predictions.csv')
hip_rcs_data <- fread('result/F1_RCS/single_right_hip_RCS_predictions.csv')
knee_rcs_data <- fread('result/F1_RCS/single_right_knee_RCS_predictions.csv')

# 排除死亡结局
double_disease <- double_rcs_data[outcome != "death"]
hip_disease    <- hip_rcs_data[outcome != "death"]
knee_disease   <- knee_rcs_data[outcome != "death"]

# 添加 trait
double_disease$trait <- "Knee - Hip"
hip_disease$trait    <- "Hip"
knee_disease$trait   <- "Knee"

# 合并数据
rcs_disease <- rbindlist(list(double_disease, hip_disease, knee_disease), use.names = T, fill = T)

# 设置 trait 顺序
rcs_disease$trait <- factor(rcs_disease$trait, levels = c("Knee", "Hip", "Knee - Hip"))

# 仅保留 3 knots
rcs_disease <- rcs_disease[knots == "3 knots"]

# 替换疾病名称
rcs_disease <- rcs_disease %>%
  mutate(
    outcome_label = dplyr::recode(
      outcome_label,
      CVD = "Cardiovascular disease",
      Respiratory = "Respiratory disease",
      Digestive = "Digestive disease",
      Endocrine = "Endocrine disease",
      Neurological = "Neurological disease",
      Psychiatric = "Psychiatric disease",
      Musculoskeletal = "Musculoskeletal disease"
    )
  )

# 设置颜色与主题
trait_colors <- c(
  "Knee" = "#006affff",
  "Hip" = "#1dcc2fff",
  "Knee - Hip" = "#FD3E5Cff"
)


my_theme <- theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 26, family = "Times New Roman", face = "bold", colour = "black"),
    axis.title = element_text(size = 24, family = "Times New Roman", face = 'bold'),
    axis.text = element_text(size = 20, family = "Times New Roman", face = 'bold', colour = 'black'),
    legend.title=element_blank(),
    legend.text = element_text(size = 18, family = "Times New Roman", face = 'bold'),
    axis.line = element_line(color  = "black", linewidth = 1.2, lineend = "square"),
    legend.position = "bottom",
    plot.margin = margin(t = 5.5, r = 2, b = 5.5, l = 2, unit = "pt")
  )


# 按疾病循环绘图
#   左：Knee；中：Hip（显示标题）；右：Knee - Hip

plot_list <- list()

disease_order <- unique(rcs_disease$outcome_label)

for (disease in disease_order) {
  
  # 提取当前疾病数据
  df_rcs <- rcs_disease[rcs_disease$outcome_label == disease,]
  
  # 当前疾病三个 trait 共用 Y 轴范围
  y_range <- range(c(df_rcs$lower, df_rcs$upper), na.rm = TRUE)
  
  y_range <- c(floor(y_range[1] * 10) / 10, ceiling(y_range[2] * 10) / 10)
  
  if (y_range[1] == y_range[2]) {
    y_range <- y_range + c(-0.1, 0.1)
  }
  
  # 左图：Knee
  df_knee <- df_rcs[df_rcs$trait == "Knee",]
  
  p1 <- ggplot(df_knee,aes(x = age_gap_adjusted, y = yhat, colour = trait, fill = trait)) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7, color = "black") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.20, colour = NA) +
    geom_line(linewidth = 1.3) +
    scale_colour_manual(values = trait_colors, breaks = names(trait_colors), drop = F) +
    scale_fill_manual(values = trait_colors,breaks = names(trait_colors), drop = F) +
    scale_x_continuous(limits = c(-10, 10),breaks = seq(-10, 10, 5)) +
    scale_y_continuous(limits = y_range) +
    labs(title = NULL,x = "Age acceleration (years)",y = "HR (95% CI)", colour = NULL, fill = NULL ) +
    my_theme
  
  # 中图：Hip
  # 显示当前疾病名称
  df_hip <- df_rcs[df_rcs$trait == "Hip",]
  
  p2 <- ggplot(df_hip, aes(x = age_gap_adjusted, y = yhat, colour = trait, fill = trait)) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7, color = "black") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.20, colour = NA) +
    geom_line(linewidth = 1.3) +
    scale_colour_manual(values = trait_colors, breaks = names(trait_colors), drop = F) +
    scale_fill_manual(values = trait_colors, breaks = names(trait_colors), drop = F) +
    scale_x_continuous(limits = c(-10, 10),breaks = seq(-10, 10, 5)) +
    scale_y_continuous(limits = y_range) +
    labs(title = disease, x = "Age acceleration (years)", y = "HR (95% CI)", colour = NULL, fill = NULL ) +
    my_theme
  
  # 右图：Knee - Hip
  df_double <- df_rcs[df_rcs$trait == "Knee - Hip",]
  
  p3 <- ggplot(df_double, aes(x = age_gap_adjusted, y = yhat, colour = trait, fill = trait)) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7, color = "black") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.20, colour = NA) +
    geom_line(linewidth = 1.3) +
    scale_colour_manual(values = trait_colors, breaks = names(trait_colors), drop = F) +
    scale_fill_manual(values = trait_colors, breaks = names(trait_colors), drop = F) +
    scale_x_continuous(limits = c(-10, 10),breaks = seq(-10, 10, 5)) +
    scale_y_continuous(limits = y_range) +
    labs(title = NULL, x = "Age acceleration (years)", y = "HR (95% CI)", colour = NULL, fill = NULL ) +
    my_theme
  
  # 拼接当前疾病的三个子图
  # 中间增加横向间距
  
  plot_list[[disease]] <- (p1 + plot_spacer() + p2 + plot_spacer() + p3) +
    plot_layout(ncol = 5, widths = c(1, 0.2, 1, 0.2, 1), guides = "collect") &
    theme(legend.position = "bottom")
}

# 4. 隐藏前面疾病行的图例，仅保留最后一行的图例

if (length(plot_list) > 1) {
  for (i in seq_len(length(plot_list) - 1)) {
    plot_list[[i]] <- plot_list[[i]] & theme(legend.position = "none",
                                             axis.title.x=element_blank())
  }
}

# 最后一行保留图例
plot_list[[length(plot_list)]] <-
  plot_list[[length(plot_list)]] & theme(legend.position = "bottom")

# 纵向拼接所有疾病
final_plot <- wrap_plots(plot_list, ncol = 1, guides = "keep")

ggsave("result/F1_RCS/rcs_right_disease.pdf", final_plot, 
       width = 16, height = 22, device = cairo_pdf, limitsize = F, bg = "white")

# RCS 左侧====
library(rms)
library(survival)
library(dplyr)
library(ggplot2)
library(data.table)
library(patchwork)

csv_files <- list.files(path = 'data', 
                        pattern = ".*_external_valid.*\\.csv$", full.names = TRUE)


for (csv_file in csv_files) {
  disease_names <- c("death", "CVD", "respiratory", "digestive", "endocrine", 
                     "neurological", "psychiatric", "musculoskeleta")
  
  disease_labels <- c(
    death = "All-cause mortality",
    CVD = "CVD",
    respiratory = "Respiratory",
    digestive = "Digestive",
    endocrine = "Endocrine",
    neurological = "Neurological",
    psychiatric = "Psychiatric",
    musculoskeleta = "Musculoskeletal"
  )
  
  
  # knots 设置
  
  knot_probs <- list(
    `3 knots` = c(0.05, 0.50, 0.95),
    `4 knots` = c(0.05, 0.35, 0.65, 0.95),
    `5 knots` = c(0.05, 0.25, 0.50, 0.75, 0.95),
    `6 knots` = c(0.05, 0.20, 0.40, 0.60, 0.80, 0.95)
  )
  
  
  # 计算 knots
  
  pred_age <- fread(csv_file)
  
  x <- pred_age$age_gap_adjusted
  
  knots_list <- lapply(knot_probs, function(p) {
    quantile(x, probs = p, na.rm = TRUE, names = FALSE)
  })
  
  knots_list
  
  # reference 和绘图范围
  
  reference_value <- median(pred_age$age_gap_adjusted, na.rm = TRUE)
  
  x_range <- quantile(pred_age$age_gap_adjusted,
                      probs = c(0.01, 0.99),
                      na.rm = TRUE)
  
  x_min <- x_range[1]
  x_max <- x_range[2]
  
  cat("Reference =", reference_value, "\n")
  cat("X-axis =", x_min, "to", x_max, "\n")
  
  
  # 输出文件夹
  
  output_dir <- "result/F1_RCS"
  
  # 设置 rms datadist
  dd <- datadist(pred_age)
  dd$limits["Adjust to", "age_gap_adjusted"] <- reference_value
  options(datadist = "dd")
  
  # 保存所有结果
  all_results <- list()
  all_predictions <- list()
  
  # 循环 knots
  for (knot_name in names(knots_list)) {
    cat("\n====================================\n")
    cat(knot_name, "\n")
    cat("====================================\n")
    
    knots <- knots_list[[knot_name]]
    
    plot_list <- list()
    
    # 10. 循环疾病
    for (disease in disease_names) {
      cat("\nRunning:", disease, "\n")
      
      # 构造 Surv
      time_var <- paste0(disease, "_time")
      event_var <- disease
      
      # 分析数据
      vars_needed <- c(
        time_var, event_var, "age_gap_adjusted", "age_i2", "sex", "Mixed", "Asian", 
        "Black", "Chinese", "others", "education_class"
      )
      
      dat <- pred_age %>%
        dplyr::select(all_of(vars_needed)) %>%
        na.omit()
      
      # 构造 formula
      fml <- as.formula(
        paste0(
          "Surv(", time_var, ", ", event_var, ") ~ ",
          "rcs(age_gap_adjusted, ", length(knots), ") + ", 
          "age_i2 + sex + Mixed + Asian + Black + Chinese + others + education_class"
        )
      )
      
      
      
      # 拟合 Cox RCS
      fit <- cph(fml, data = dat, x = TRUE, y = TRUE, surv = TRUE)
      
      # ANOVA
      a <- anova(fit)
      # AIC
      model_aic <- AIC(fit)
      # 提取 P overall / P nonlinear
      a_df <- as.data.frame(a)
      row_names <- rownames(a_df)
      
      # 找到 age_gap_adjusted 对应行
      overall_row <- which(grepl("^age_gap_adjusted$", row_names))
      nonlinear_row <- which(grepl("Nonlinear", row_names, ignore.case = TRUE))
      
      # 如果找不到，设置 NA
      p_overall <- NA_real_
      p_nonlinear <- NA_real_
      
      if (length(overall_row) > 0) {
        p_col <- grep("P|Wald", colnames(a_df), ignore.case = TRUE, value = TRUE)
        if (length(p_col) > 0) {
          p_overall <- a_df[overall_row[1], p_col[length(p_col)]]
        }
      }
      
      if (length(nonlinear_row) > 0) {
        p_col <- grep("P|Wald", colnames(a_df), ignore.case = TRUE, value = TRUE)
        if (length(p_col) > 0) {
          p_nonlinear <- a_df[nonlinear_row[1], p_col[length(p_col)]]
        }
      }
      
      # Prediction
      pred <- rms::Predict(
        fit, age_gap_adjusted = seq(x_min, x_max, length.out = 300),
        ref.zero = TRUE, fun = exp, conf.int = 0.95
      )
      pred_df <- as.data.frame(pred)
      
      # 保存 prediction
      pred_df <- pred_df %>%
        mutate(
          outcome = disease,
          outcome_label = disease_labels[disease],
          knots = knot_name,
          reference = reference_value,
          N = nrow(dat),
          events = sum(dat[[event_var]] == 1, na.rm = TRUE),
          P_overall = p_overall,
          P_nonlinear = p_nonlinear
        )
      
      all_predictions[[paste0(knot_name, "_", disease)]] <- pred_df
      
      # P value 格式
      p_overall_text <- ifelse(
        is.na(p_overall),
        "P overall = NA",
        ifelse(
          p_overall < 0.001,
          "P overall < 0.001",
          paste0("P overall = ", formatC(
            p_overall, format = "f", digits = 3
          ))
        )
      )
      
      p_nonlinear_text <- ifelse(
        is.na(p_nonlinear),
        "P nonlinear = NA",
        ifelse(
          p_nonlinear < 0.001,
          "P nonlinear < 0.001",
          paste0("P nonlinear = ", formatC(
            p_nonlinear, format = "f", digits = 3
          ))
        )
      )
      
      # 单个 disease 图
      p <- ggplot(pred_df, aes(x = age_gap_adjusted, y = yhat)) +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.20) +
        geom_line(linewidth = 0.9) +
        geom_hline(yintercept = 1,
                   linetype = "dashed",
                   linewidth = 0.5) +
        geom_vline(xintercept = reference_value,
                   linetype = "dotted",
                   linewidth = 0.5) +
        labs(
          title = disease_labels[disease],
          subtitle = paste(p_overall_text, "\n", p_nonlinear_text),
          x = "Age gap (adjusted)",
          y = "Hazard ratio (95% CI)"
        ) +
        coord_cartesian(xlim = c(x_min, x_max)) +
        theme_classic(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9),
          axis.title = element_text(face = "bold")
        )
      plot_list[[disease]] <- p
      
      # 保存统计结果
      all_results[[paste0(knot_name, "_", disease)]] <-
        data.frame(
          outcome = disease,
          outcome_label = disease_labels[disease],
          knots = knot_name,
          aic = model_aic,
          knot_number = length(knots),
          knot_positions = paste(round(knots, 4), collapse = ", "),
          N = nrow(dat),
          events = sum(dat[[event_var]] == 1, na.rm = TRUE),
          reference = reference_value,
          x_min = x_min,
          x_max = x_max,
          P_overall = p_overall,
          P_nonlinear = p_nonlinear
        )
    }
    
    # 8个疾病合并成一张图
    combined_plot <- wrap_plots(plot_list, ncol = 2) +
      plot_annotation(
        title = paste0("Restricted cubic spline analysis: ", knot_name),
        subtitle = paste0("Reference = median age gap (", round(reference_value, 3), ")")
      )
    
    # 保存图片
    file_name <- paste0(gsub("external_valid.csv", "", basename(csv_file)), gsub(" ", "_", knot_name))
    
    ggsave(
      filename = file.path(output_dir, paste0("RCS_", file_name, ".png")),
      plot = combined_plot, width = 7, height = 14, dpi = 300
    )
    
    # 同时保存 PDF
    ggsave(
      filename = file.path(output_dir, paste0("RCS_", file_name, ".pdf")),
      plot = combined_plot, width = 7, height = 14
    )
    cat("\nSaved:", file.path(output_dir, paste0("RCS_", file_name, ".png")), "\n")
  }
  
  # 13. 合并统计结果
  results_df <- bind_rows(all_results)
  results_df <- results_df %>%
    arrange(knot_number, match(outcome, disease_names))
  
  
  # 保存 summary
  write.csv(results_df,
            file.path(output_dir, file_name <- paste0(gsub("external_valid.csv", "", basename(csv_file)), "RCS_summary.csv")),
            row.names = FALSE)
  
  
  # 合并所有 prediction
  prediction_df <- bind_rows(all_predictions)
  write.csv(prediction_df,
            file.path(output_dir, file_name <- paste0(gsub("external_valid.csv", "", basename(csv_file)), "RCS_predictions.csv")),
            row.names = FALSE)
  
  # 查看结果
  print(results_df)
}

## 绘图 death====
double_rcs_data <- fread('result/F1_RCS/double_left_RCS_predictions.csv')
hip_rcs_data <- fread('result/F1_RCS/single_left_hip_RCS_predictions.csv')
knee_rcs_data <- fread('result/F1_RCS/single_left_knee_RCS_predictions.csv')
head(double_rcs_data)

double_death <- double_rcs_data[outcome == "death"]
hip_death    <- hip_rcs_data[outcome == "death"]
knee_death   <- knee_rcs_data[outcome == "death"]

# 添加 trait
double_death$trait <- 'Knee - Hip'
hip_death$trait <- 'Hip'
knee_death$trait <- 'Knee'

rcs_death <- rbindlist(list(double_death, hip_death, knee_death),
                       use.names = TRUE,
                       fill = TRUE)

rcs_death$trait <- factor(rcs_death$trait, levels = c("Knee", "Hip", "Knee - Hip"))

rcs_death <- rcs_death[rcs_death$knots == '3 knots',]

# 统一坐标轴范围
y_range <- range(c(rcs_death$lower, rcs_death$upper), na.rm = TRUE)

# 如果希望上下留一点空白
y_range <- c(floor(y_range[1] * 10) / 10, ceiling(y_range[2] * 10) / 10)

# 统一主题
my_theme <- theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 26, family = "Times New Roman", face = "bold", colour = "black"),
    axis.title = element_text(size = 24, family = "Times New Roman", face = 'bold'),
    axis.text = element_text(size = 20, family = "Times New Roman", face = 'bold', colour = 'black'),
    legend.title=element_blank(),
    legend.text = element_text(size = 18, family = "Times New Roman", face = 'bold'),
    axis.line = element_line(color  = "black", linewidth = 1.2, lineend = "square")
  )


# Knee
d1 <- rcs_death[rcs_death$trait == "Knee", ]

p1 <- ggplot(d1, aes(x = age_gap_adjusted, y = yhat)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7, color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#006affff", alpha = 0.20) +
  geom_line(color = "#006affff", linewidth = 1.3) +
  scale_x_continuous(limits = c(-9,9), breaks = seq(-9,9,3)) +
  scale_y_continuous(limits = y_range) +
  labs(title = "Knee", x = "Age acceleration (years)", y = "Hazard ratio (95% CI)") +
  my_theme

p1

# Hip
d2 <- rcs_death[rcs_death$trait == "Hip", ]

p2 <- ggplot(d2, aes(x = age_gap_adjusted, y = yhat)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7, color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#1dcc2fff", alpha = 0.20) +
  geom_line(color = "#1dcc2fff", linewidth = 1.3) +
  scale_x_continuous(limits = c(-9,9), breaks = seq(-9,9,3)) +
  scale_y_continuous(limits = y_range) +
  labs(title = "Hip", x = "Age acceleration (years)", y = "Hazard ratio (95% CI)") +
  my_theme

# Knee - Hip
d3 <- rcs_death[rcs_death$trait == "Knee - Hip", ]

p3 <- ggplot(d3, aes(x = age_gap_adjusted, y = yhat)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7, color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#FD3E5Cff", alpha = 0.20) +
  geom_line(color = "#FD3E5Cff", linewidth = 1.3) +
  scale_x_continuous(limits = c(-9,9), breaks = seq(-9,9,3)) +
  scale_y_continuous(limits = y_range) +
  labs(title = "Knee - Hip", x = "Age acceleration (years)", y = "Hazard ratio (95% CI)") +
  my_theme



p <- p1 + p2 + p3 + plot_layout(ncol = 3)
p

ggsave('result/F1_RCS/rcs_left_knee.png', p1, width = 6, height = 6,  dpi = 600)
ggsave('result/F1_RCS/rcs_left_hip.png', p2, width = 6, height = 6,  dpi = 600)
ggsave('result/F1_RCS/rcs_left_double.png', p3, width = 6, height = 6,  dpi = 600)
ggsave('result/F1_RCS/rcs_left_all.png', p, width = 18, height = 6,  dpi = 600)

## 绘图 disease ====

double_rcs_data <- fread('result/F1_RCS/double_left_RCS_predictions.csv')
hip_rcs_data <- fread('result/F1_RCS/single_left_hip_RCS_predictions.csv')
knee_rcs_data <- fread('result/F1_RCS/single_left_knee_RCS_predictions.csv')

# 排除死亡结局
double_disease <- double_rcs_data[outcome != "death"]
hip_disease    <- hip_rcs_data[outcome != "death"]
knee_disease   <- knee_rcs_data[outcome != "death"]

# 添加 trait
double_disease$trait <- "Knee - Hip"
hip_disease$trait    <- "Hip"
knee_disease$trait   <- "Knee"

# 合并数据
rcs_disease <- rbindlist(list(double_disease, hip_disease, knee_disease), use.names = T, fill = T)

# 设置 trait 顺序
rcs_disease$trait <- factor(rcs_disease$trait, levels = c("Knee", "Hip", "Knee - Hip"))

# 仅保留 3 knots
rcs_disease <- rcs_disease[knots == "3 knots"]

# 替换疾病名称
rcs_disease <- rcs_disease %>%
  mutate(
    outcome_label = dplyr::recode(
      outcome_label,
      CVD = "Cardiovascular disease",
      Respiratory = "Respiratory disease",
      Digestive = "Digestive disease",
      Endocrine = "Endocrine disease",
      Neurological = "Neurological disease",
      Psychiatric = "Psychiatric disease",
      Musculoskeletal = "Musculoskeletal disease"
    )
  )

# 设置颜色与主题
trait_colors <- c(
  "Knee" = "#006affff",
  "Hip" = "#1dcc2fff",
  "Knee - Hip" = "#FD3E5Cff"
)


my_theme <- theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 26, family = "Times New Roman", face = "bold", colour = "black"),
    axis.title = element_text(size = 24, family = "Times New Roman", face = 'bold'),
    axis.text = element_text(size = 20, family = "Times New Roman", face = 'bold', colour = 'black'),
    legend.title=element_blank(),
    legend.text = element_text(size = 18, family = "Times New Roman", face = 'bold'),
    axis.line = element_line(color  = "black", linewidth = 1.2, lineend = "square"),
    legend.position = "bottom",
    plot.margin = margin(t = 5.5, r = 2, b = 5.5, l = 2, unit = "pt")
  )


# 按疾病循环绘图
#   左：Knee；中：Hip（显示标题）；右：Knee - Hip

plot_list <- list()

disease_order <- unique(rcs_disease$outcome_label)

for (disease in disease_order) {
  
  # 提取当前疾病数据
  df_rcs <- rcs_disease[rcs_disease$outcome_label == disease,]
  
  # 当前疾病三个 trait 共用 Y 轴范围
  y_range <- range(c(df_rcs$lower, df_rcs$upper), na.rm = TRUE)
  
  y_range <- c(floor(y_range[1] * 10) / 10, ceiling(y_range[2] * 10) / 10)
  
  if (y_range[1] == y_range[2]) {
    y_range <- y_range + c(-0.1, 0.1)
  }
  
  # 左图：Knee
  df_knee <- df_rcs[df_rcs$trait == "Knee",]
  
  p1 <- ggplot(df_knee,aes(x = age_gap_adjusted, y = yhat, colour = trait, fill = trait)) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7, color = "black") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.20, colour = NA) +
    geom_line(linewidth = 1.3) +
    scale_colour_manual(values = trait_colors, breaks = names(trait_colors), drop = F) +
    scale_fill_manual(values = trait_colors,breaks = names(trait_colors), drop = F) +
    scale_x_continuous(limits = c(-10, 10),breaks = seq(-10, 10, 5)) +
    scale_y_continuous(limits = y_range) +
    labs(title = NULL,x = "Age acceleration (years)",y = "HR (95% CI)", colour = NULL, fill = NULL ) +
    my_theme
  
  # 中图：Hip
  # 显示当前疾病名称
  df_hip <- df_rcs[df_rcs$trait == "Hip",]
  
  p2 <- ggplot(df_hip, aes(x = age_gap_adjusted, y = yhat, colour = trait, fill = trait)) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7, color = "black") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.20, colour = NA) +
    geom_line(linewidth = 1.3) +
    scale_colour_manual(values = trait_colors, breaks = names(trait_colors), drop = F) +
    scale_fill_manual(values = trait_colors, breaks = names(trait_colors), drop = F) +
    scale_x_continuous(limits = c(-10, 10),breaks = seq(-10, 10, 5)) +
    scale_y_continuous(limits = y_range) +
    labs(title = disease, x = "Age acceleration (years)", y = "HR (95% CI)", colour = NULL, fill = NULL ) +
    my_theme
  
  # 右图：Knee - Hip
  df_double <- df_rcs[df_rcs$trait == "Knee - Hip",]
  
  p3 <- ggplot(df_double, aes(x = age_gap_adjusted, y = yhat, colour = trait, fill = trait)) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7, color = "black") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.20, colour = NA) +
    geom_line(linewidth = 1.3) +
    scale_colour_manual(values = trait_colors, breaks = names(trait_colors), drop = F) +
    scale_fill_manual(values = trait_colors, breaks = names(trait_colors), drop = F) +
    scale_x_continuous(limits = c(-10, 10),breaks = seq(-10, 10, 5)) +
    scale_y_continuous(limits = y_range) +
    labs(title = NULL, x = "Age acceleration (years)", y = "HR (95% CI)", colour = NULL, fill = NULL ) +
    my_theme
  
  # 拼接当前疾病的三个子图
  # 中间增加横向间距
  
  plot_list[[disease]] <- (p1 + plot_spacer() + p2 + plot_spacer() + p3) +
    plot_layout(ncol = 5, widths = c(1, 0.2, 1, 0.2, 1), guides = "collect") &
    theme(legend.position = "bottom")
}

# 4. 隐藏前面疾病行的图例，仅保留最后一行的图例

if (length(plot_list) > 1) {
  for (i in seq_len(length(plot_list) - 1)) {
    plot_list[[i]] <- plot_list[[i]] & theme(legend.position = "none",
                                             axis.title.x=element_blank())
  }
}

# 最后一行保留图例
plot_list[[length(plot_list)]] <-
  plot_list[[length(plot_list)]] & theme(legend.position = "bottom")

# 纵向拼接所有疾病
final_plot <- wrap_plots(plot_list, ncol = 1, guides = "keep")

ggsave("result/F1_RCS/rcs_left_disease.pdf", final_plot, 
       width = 16, height = 22, device = cairo_pdf, limitsize = F, bg = "white")

# 预期寿命数据====
csv_files <- list.files(path = 'data', 
                        pattern = ".*oof.*\\.csv$", full.names = TRUE)
csv_files <- list.files(path = 'data', 
                        pattern = ".*(oof|valid).*\\.csv$", 
                        full.names = TRUE,
                        ignore.case = TRUE)

for (csv_file in csv_files) {
  data <- fread(csv_file)
  data <- data %>%
    mutate(
      death_age = age_i2 + death_time / 365.25,
      decile4 = ntile(age_gap_adjusted, 4),
      decile10 = ntile(age_gap_adjusted, 10)
    )
  death_data <- data %>%
    dplyr::select(
      eid, death, death_time, death_age, age_i2, age_gap_adjusted, 
      sex, education_class, Black, Asian, Mixed, others, Chinese,
      decile4, decile10
    )
  print(summary(death_data$age_i2))
  summary(death_data$death_time)
  summary(death_data$death_age)
  
  death_data <- na.omit(death_data)
  
  file_name <- gsub('oof', 'death', csv_file)
  file_name <- gsub('external_valid', 'death', file_name)
  print(file_name)
  write.csv(death_data, file_name,row.names = F)
}

# 绘图
library(haven)

dta_files <- list.files(path = 'result/F1_预期寿命', 
                        pattern = "*4.dta$", full.names = TRUE)

dta_files

## left ====
left_data <- map_dfr(dta_files[grepl("left", dta_files)], ~ read_dta(.x) %>%
                       mutate(
                         trait = case_when(
                           grepl("double_left", .x) ~ "Knee - Hip",
                           grepl("single_left_hip", .x) ~ "Hip",
                           grepl("single_left_knee", .x) ~ "Knee"))  %>%
                       mutate(trait = factor(trait, levels = c("Knee", "Hip", "Knee - Hip")))
)


# 统一坐标轴范围
y_range <- range(c(left_data$erl_lci, left_data$erl_uci), na.rm = TRUE)

# 如果希望上下留一点空白
y_range <- c(floor(y_range[1] * 10) / 10, ceiling(y_range[2] * 10) / 10)

# 统一主题
my_theme <- theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 26, family = "Times New Roman", face = "bold", colour = "black"),
    axis.title = element_text(size = 24, family = "Times New Roman", face = 'bold'),
    axis.text = element_text(size = 20, family = "Times New Roman", face = 'bold', colour = 'black'),
    legend.title=element_blank(),
    legend.text = element_text(size = 18, family = "Times New Roman", face = 'bold'),
    axis.line = element_line(color  = "black", linewidth = 1.2, lineend = "square")
  )

# Knee
d1 <- left_data[left_data$trait == "Knee", ]

p1 <- ggplot(d1, aes(x = cond, y = erl_mm)) +
  geom_ribbon(aes(ymin = erl_lci, ymax = erl_uci),  fill = "#006affff", alpha = 0.15) +   
  geom_line(color = "#006affff", linewidth = 1.2) +       
  labs(title = 'Knee', x = "Age (years)", y = "Estimated years of life lost") +
  scale_x_continuous(limits = c(45, 100), breaks = seq(45, 100, 5)) +
  scale_y_continuous(limits = y_range) +
  my_theme

p1


# Hip
d2 <- left_data[left_data$trait == "Hip", ]

p2 <- ggplot(d2, aes(x = cond, y = erl_mm)) +
  geom_ribbon(aes(ymin = erl_lci, ymax = erl_uci),  fill = "#1dcc2fff", alpha = 0.15) +   
  geom_line(color = "#1dcc2fff", linewidth = 1.2) +       
  labs(title = 'Hip', x = "Age (years)", y = "Estimated years of life lost") +
  scale_x_continuous(limits = c(45, 100), breaks = seq(45, 100, 5)) +
  scale_y_continuous(limits = y_range) +
  my_theme

p2


# Knee - Hip
d3 <- left_data[left_data$trait == "Knee - Hip", ]

p3 <- ggplot(d3, aes(x = cond, y = erl_mm)) +
  geom_ribbon(aes(ymin = erl_lci, ymax = erl_uci),  fill = "#FD3E5Cff", alpha = 0.15) +   
  geom_line(color = "#FD3E5Cff", linewidth = 1.2) +       
  labs(title = 'Knee - Hip', x = "Age (years)", y = "Estimated years of life lost") +
  scale_x_continuous(limits = c(45, 100), breaks = seq(45, 100, 5)) +
  scale_y_continuous(limits = y_range) +
  my_theme

p3


p <- p1 + p2 + p3 + plot_layout(ncol = 3)
p

ggsave('result/F1_预期寿命/life_left_knee.png', p1, width = 6, height = 6,  dpi = 600)
ggsave('result/F1_预期寿命/life_left_hip.png', p2, width = 6, height = 6,  dpi = 600)
ggsave('result/F1_预期寿命/life_left_double.png', p3, width = 6, height = 6,  dpi = 600)
ggsave('result/F1_预期寿命/life_left_all.png', p, width = 18, height = 6,  dpi = 600)



## right ====
right_data <- map_dfr(dta_files[grepl("right", dta_files)], ~ read_dta(.x) %>%
                        mutate(
                          trait = case_when(
                            grepl("double_right", .x) ~ "Knee - Hip",
                            grepl("single_right_hip", .x) ~ "Hip",
                            grepl("single_right_knee", .x) ~ "Knee"))  %>%
                        mutate(trait = factor(trait, levels = c("Knee", "Hip", "Knee - Hip")))
)


# 统一坐标轴范围
y_range <- range(c(right_data$erl_lci, right_data$erl_uci), na.rm = TRUE)

# 如果希望上下留一点空白
y_range <- c(floor(y_range[1] * 10) / 10, ceiling(y_range[2] * 10) / 10)

# 统一主题
my_theme <- theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 26, family = "Times New Roman", face = "bold", colour = "black"),
    axis.title = element_text(size = 24, family = "Times New Roman", face = 'bold'),
    axis.text = element_text(size = 20, family = "Times New Roman", face = 'bold', colour = 'black'),
    legend.title=element_blank(),
    legend.text = element_text(size = 18, family = "Times New Roman", face = 'bold'),
    axis.line = element_line(color  = "black", linewidth = 1.2, lineend = "square")
  )

# Knee
d1 <- right_data[right_data$trait == "Knee", ]

p1 <- ggplot(d1, aes(x = cond, y = erl_mm)) +
  geom_ribbon(aes(ymin = erl_lci, ymax = erl_uci),  fill = "#006affff", alpha = 0.15) +   
  geom_line(color = "#006affff", linewidth = 1.2) +       
  labs(title = "Knee", x = "Age (years)", y = "Estimated years of life lost") +
  scale_x_continuous(limits = c(45, 100), breaks = seq(45, 100, 5)) +
  scale_y_continuous(limits = y_range) +
  my_theme

p1


# Hip
d2 <- right_data[right_data$trait == "Hip", ]

p2 <- ggplot(d2, aes(x = cond, y = erl_mm)) +
  geom_ribbon(aes(ymin = erl_lci, ymax = erl_uci),  fill = "#1dcc2fff", alpha = 0.15) +   
  geom_line(color = "#1dcc2fff", linewidth = 1.2) +       
  labs(title = "Hip", x = "Age (years)", y = "Estimated years of life lost") +
  scale_x_continuous(limits = c(45, 100), breaks = seq(45, 100, 5)) +
  scale_y_continuous(limits = y_range) +
  my_theme

p2


# Knee - Hip
d3 <- right_data[right_data$trait == "Knee - Hip", ]

p3 <- ggplot(d3, aes(x = cond, y = erl_mm)) +
  geom_ribbon(aes(ymin = erl_lci, ymax = erl_uci),  fill = "#FD3E5Cff", alpha = 0.15) +   
  geom_line(color = "#FD3E5Cff", linewidth = 1.2) +       
  labs(title = "Knee - Hip", x = "Age (years)", y = "Estimated years of life lost") +
  scale_x_continuous(limits = c(45, 100), breaks = seq(45, 100, 5)) +
  scale_y_continuous(limits = y_range) +
  my_theme

p3


p <- p1 + p2 + p3 + plot_layout(ncol = 3)
p

ggsave('result/F1_预期寿命/life_right_knee.png', p1, width = 6, height = 6,  dpi = 600)
ggsave('result/F1_预期寿命/life_right_hip.png', p2, width = 6, height = 6,  dpi = 600)
ggsave('result/F1_预期寿命/life_right_double.png', p3, width = 6, height = 6,  dpi = 600)
ggsave('result/F1_预期寿命/life_right_all.png', p, width = 18, height = 6,  dpi = 600)




# 偏相关====
library(ppcor)
# 指标
cor <- fread('data/cor_data_i0.csv')
names(cor)
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
csv_files <- list.files(path = 'data', 
                        pattern = ".*(oof|valid).*\\.csv$", 
                        full.names = TRUE,
                        ignore.case = TRUE)
for (csv_file in csv_files) {
  pred_age <- fread(csv_file)
  
  data <- pred_age %>%
    inner_join(cor, by = "eid") %>%
    inner_join(p53[,c('eid','followup_time')], by = "eid") %>%
    inner_join(age_v0, by = "eid") 
  names(data)
  
  # follow-up time 转成年
  data <- data %>%mutate(followup_year = followup_time / 365.25)
  
  # V0 指标
  markers <- c('father_death_age','mother_death_age','tl','SII',
               'CRP','frailty_score','BMD','grip_min','grip_max')
  
  # 偏相关校正变量
  covariates <- c("age_i2", "followup_year", "sex", "Mixed", "Asian", "Black", "Chinese", "others", "education_class")
  
  # 批量偏相关
  results <- map_dfr(markers, function(marker) {
    tmp <- data %>%
      dplyr::select(all_of(marker), age_gap_adjusted, all_of(covariates)) %>%
      na.omit() %>%
      as.data.frame()
    
    pc <- pcor.test(x = tmp[[marker]], y = tmp$age_gap_adjusted, z = tmp[, covariates], method = "pearson")
    
    tibble(marker = marker, N = nrow(tmp), partial_r = as.numeric(pc$estimate), 
           partial_r2 = as.numeric(pc$estimate^2),p = as.numeric(pc$p.value))
  })
  file_name <- paste0('result/F1_偏相关/acc_', gsub('_oof', '', basename(csv_file)))
  write.csv(results, file_name,row.names = F)
  
  # 偏相关校正变量--去除年龄
  covariates <- c("followup_year", "sex", "Mixed", "Asian", "Black", "Chinese", "others", "education_class")
  
  # 批量偏相关
  results <- map_dfr(markers, function(marker) {
    tmp <- data %>%
      dplyr::select(all_of(marker), adjusted_age, all_of(covariates)) %>%
      na.omit() %>%
      as.data.frame()
    
    pc <- pcor.test(x = tmp[[marker]], y = tmp$adjusted_age, z = tmp[, covariates], method = "pearson")
    
    tibble(marker = marker, N = nrow(tmp), partial_r = as.numeric(pc$estimate), 
           partial_r2 = as.numeric(pc$estimate^2),p = as.numeric(pc$p.value))
  })
  
  file_name <- paste0('result/F1_偏相关/age_', gsub('_oof', '', basename(csv_file)))
  file_name <- gsub('_external_valid', '', file_name)
  write.csv(results, file_name, row.names = F)
  
}


library(tidyverse)

# 文件路径
path <- "/path/to/epi/result/F1_偏相关"
files <- c(
  "acc_single_right_knee.csv", "acc_single_right_hip.csv",
  "acc_double_right.csv", "age_single_right_knee.csv",
  "age_single_right_hip.csv", "age_double_right.csv",
  "acc_single_left_knee.csv", "acc_single_left_hip.csv",
  "acc_double_left.csv", "age_single_left_knee.csv",
  "age_single_left_hip.csv", "age_double_left.csv"
)

# 2. 批量读取
all_results <- map_dfr(files, function(file) {
  df <- read.csv(file.path(path, file), stringsAsFactors = FALSE)
  
  outcome <- ifelse(grepl("^acc_", file), "Age acceleration", "Predicted age")
  model <- case_when(
    grepl("single_right_knee", file) ~ "Single right knee",
    grepl("single_right_hip", file)  ~ "Single right hip",
    grepl("double_right", file)      ~ "Double right",
    grepl("single_left_knee", file) ~ "Single left knee",
    grepl("single_left_hip", file)  ~ "Single left hip",
    grepl("double_left", file)      ~ "Double left",
    TRUE ~ NA_character_
  )
  
  df %>%
    mutate(outcome = outcome, model = model)
})


# 3. 显著性标记
all_results <- all_results %>%
  filter(!marker %in% c('CRP', 'SII', 'grip_min')) %>%
  mutate(significance = case_when(p < 0.001 ~ "***", p < 0.01  ~ "**", p < 0.05  ~ "*", TRUE      ~ ""))

# 4. 指标顺序
marker_mapping <- c(
  "father_death_age" = "Father's age at death",
  "mother_death_age" = "Mother's age at death",
  "tl" = "Telomere length",
  "frailty_score" = "Frailty score",
  "BMD" = "Heel bone nineral density",
  "grip_max" = "Hand grip strength"
)

model_mapping <- c(
  "Single right knee" = "Knee",
  "Single right hip" = "Hip",
  "Double right" = "Knee - Hip",
  "Single left knee" = "Knee",
  "Single left hip" = "Hip",
  "Double left" = "Knee - Hip"
)

all_results <- all_results %>%
  mutate(marker = factor(marker, levels = names(marker_mapping), labels = unname(marker_mapping)),
         trait = factor(model, levels = names(model_mapping), labels = unname(model_mapping))) %>%
  mutate(marker = factor(marker,
                         levels =  rev(c("Father's age at death", "Mother's age at death", 
                                         "Telomere length", 'Heel bone nineral density', 
                                         "Hand grip strength","Frailty score"))), 
         trait = factor(trait, levels = c("Knee", "Hip", "Knee - Hip"))) 




library(tidyverse)

# 3. Age acceleration
acc_data <- all_results %>%
  filter(outcome == "Age acceleration") %>%
  filter(grepl('right', model))

p_acc <- ggplot(acc_data, aes(x = trait, y = marker, fill = partial_r)) +
  geom_tile(color = "white", linewidth = 0.5,width = 0.98, height = 0.98) +
  geom_text(aes(label = significance), size = 9) +
  scale_fill_gradient2(
    low = "#108b96", mid = "white", high = "#ff7e0f",
    midpoint = 0, limits = c(-0.1, 0.1),
    name = "Correlation coefficients (r)",
    guide = guide_colorbar(
      title.position = "left",        # 标题放到色条左侧
      title.theme = element_text(
        angle = 90, hjust = 0.5,  vjust = 0.5, family = "Times New Roman",
        size = 18, face = "bold", colour = "black"
        ),
      barwidth  = unit(0.8, "cm"),
      barheight = unit(10, "cm")      # ← 色条加长
    ))+
  labs(title = NULL, x = NULL, y = NULL) +
  theme_minimal() + 
  theme(
    plot.title = element_text(family='Times New Roman', size = 26,  face = 'bold', vjust = 1),
    legend.title = element_text(size = 20, family = "Times New Roman", colour = 'black', face = "bold"),
    legend.text = element_text(family='Times New Roman', size = 16, face = 'bold'),
    axis.title.y = element_blank(), # 不显示Y轴标题
    axis.title.x = element_blank(), # 不显示X轴标题
    panel.grid = element_blank(),  # 不显示网格线
    plot.background = element_rect(fill = "white", color = "white"),  # 设置整体背景为白色
    panel.background = element_rect(fill = "white", color = "white"),  # 设置面板背景为白色
    axis.text.x = element_text(size = 22, family = "Times New Roman", colour = 'black', face = "bold"),
    axis.text.y = element_text(size = 22, family = "Times New Roman", colour = 'black', face = "bold"),
    legend.key.size = unit(1, "cm")
  ) +
  coord_fixed(ratio = 1)
p_acc

ggsave('result/F1_偏相关/right_acc.pdf', p_acc, height = 12, width = 9, device = cairo_pdf)
ggsave('result/F1_偏相关/right_acc.png', p_acc, height = 12, width = 9, dpi = 600)
# 4. Predicted age
age_data <- all_results %>%
  filter(outcome == "Predicted age")  %>%
  filter(marker != 'CRP')


p_age <- ggplot(age_data, aes(x = model, y = marker, fill = partial_r)) +
  geom_tile(color = "white", linewidth = 0.5,width = 0.98, height = 0.98) +
  geom_text(aes(label = significance), size = 9) +
  scale_fill_gradient2(low = "#108b96",mid = "white",high = "#ff7e0f",
                       midpoint = 0,  limits = c(-0.26, 0.26), 
                       breaks = c(-0.26,-0.13,0,0.13, 0.26),
                       name = "Correlation coefficients (r)",
                       guide = guide_colorbar(
                         title.position = "left",        # 标题放到色条左侧
                         title.theme = element_text(
                           angle = 90, hjust = 0.5,  vjust = 0.5, family = "Times New Roman",
                           size = 18, face = "bold", colour = "black"
                         ),
                         barwidth  = unit(0.8, "cm"),
                         barheight = unit(10, "cm")      # ← 色条加长
                       ))+
  labs(title = NULL, x = NULL, y = NULL) +
  theme_minimal() + 
  theme(
    plot.title = element_text(family='Times New Roman', size = 26,  face = 'bold', vjust = 1),
    legend.title = element_text(size = 20, family = "Times New Roman", colour = 'black', face = "bold"),
    legend.text = element_text(family='Times New Roman', size = 16, face = 'bold'),
    axis.title.y = element_blank(), # 不显示Y轴标题
    axis.title.x = element_blank(), # 不显示X轴标题
    panel.grid = element_blank(),  # 不显示网格线
    plot.background = element_rect(fill = "white", color = "white"),  # 设置整体背景为白色
    panel.background = element_rect(fill = "white", color = "white"),  # 设置面板背景为白色
    axis.text.x = element_text(size = 22, family = "Times New Roman", colour = 'black', face = "bold"),
    axis.text.y = element_text(size = 22, family = "Times New Roman", colour = 'black', face = "bold"),
    legend.key.size = unit(1, "cm")
  ) +
  coord_fixed(ratio = 1)

# 查看
p_age

ggsave('result/F1_偏相关/right_age.pdf', p_age, height = 15, width = 9, device = cairo_pdf)
ggsave('result/F1_偏相关/right_age.png', p_age, height = 15, width = 9, dpi = 600)


# 5. Left age acceleration
acc_data <- all_results %>%
  filter(outcome == "Age acceleration") %>%
  filter(grepl('left', model))

p_acc <- ggplot(acc_data, aes(x = trait, y = marker, fill = partial_r)) +
  geom_tile(color = "white", linewidth = 0.5,width = 0.98, height = 0.98) +
  geom_text(aes(label = significance), size = 9) +
  scale_fill_gradient2(
    low = "#108b96", mid = "white", high = "#ff7e0f",
    midpoint = 0, limits = c(-0.1, 0.1),
    name = "Correlation coefficients (r)",
    guide = guide_colorbar(
      title.position = "left",        # 标题放到色条左侧
      title.theme = element_text(
        angle = 90, hjust = 0.5,  vjust = 0.5, family = "Times New Roman",
        size = 18, face = "bold", colour = "black"
      ),
      barwidth  = unit(0.8, "cm"),
      barheight = unit(10, "cm")      # ← 色条加长
    ))+
  labs(title = NULL, x = NULL, y = NULL) +
  theme_minimal() + 
  theme(
    plot.title = element_text(family='Times New Roman', size = 26,  face = 'bold', vjust = 1),
    legend.title = element_text(size = 20, family = "Times New Roman", colour = 'black', face = "bold"),
    legend.text = element_text(family='Times New Roman', size = 16, face = 'bold'),
    axis.title.y = element_blank(), # 不显示Y轴标题
    axis.title.x = element_blank(), # 不显示X轴标题
    panel.grid = element_blank(),  # 不显示网格线
    plot.background = element_rect(fill = "white", color = "white"),  # 设置整体背景为白色
    panel.background = element_rect(fill = "white", color = "white"),  # 设置面板背景为白色
    axis.text.x = element_text(size = 22, family = "Times New Roman", colour = 'black', face = "bold"),
    axis.text.y = element_text(size = 22, family = "Times New Roman", colour = 'black', face = "bold"),
    legend.key.size = unit(1, "cm")
  ) +
  coord_fixed(ratio = 1)
p_acc

ggsave('result/F1_偏相关/left_acc.pdf', p_acc, height = 12, width = 9, device = cairo_pdf)
ggsave('result/F1_偏相关/left_acc.png', p_acc, height = 12, width = 9, dpi = 600)

