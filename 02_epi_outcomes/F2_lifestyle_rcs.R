# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 02_epi_outcomes
# 生活方式因素与年龄加速的限制性立方样条（RCS）分析
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
# 生活方式 RCS====
library(data.table)
library(dplyr)
library(rms)
library(ggplot2)
library(patchwork)

# 随访时间 V0-V2
p53 <- fread('shared_project/data/p53.csv')
p53 <- p53 %>%
  mutate(across(where(is.character), ~ if_else(str_detect(., "^\\s*$"), NA_character_, .))) %>%
  mutate(followup_time = as.numeric(as.Date(p53_i2) - as.Date(p53_i0)))

# V0 age
age_v0 <- fread('shared_project/data/sex_age_race.csv')
age_v0 <- age_v0[,c(1,3)]

# 生活方式
life_data <- fread('data/life.csv')
lifestyle_vars <- c(
  "sleep_duration", "Overall_acceleration_average_valid",
  "MET_minutes_per_week_for_moderate_activity", "MET_minutes_per_week_for_vigorous_activity",
  "MET_minutes_per_week_for_walking", "Summed_MET_minutes_per_week_for_all_activity",
  "bmi", "smoke_pack", "somke_prop", "drink", "diet_i0_score", "AMED_score", "AHEI_score"
)

table(life_data$sleep_duration)
life_data$sleep_duration <- as.numeric(life_data$sleep_duration)
# 年龄加速
files <- c(
  double_right = "data/double_right_oof.csv",
  single_right_hip = "data/single_right_hip_oof.csv",
  single_right_knee = "data/single_right_knee_oof.csv"
)

# 协变量
covariates <- c("age_i0", "followup_year", "sex", "Mixed", "Asian", "Black", 
                "Chinese", "others", "education_class")

# 结果
AIC_result <- list()
pred_result <- list()
# 输出文件夹
out_dir <- "result/F2_生活方式RCS"

# 循环年龄加速
for (outcome in names(files)) {
  message("\n================ ", outcome, " ================\n")
  
  d <- fread(files[outcome]) %>%
    inner_join(life_data, by = "eid") %>%
    inner_join(p53[, c("eid", "followup_time")], by = "eid") %>%
    inner_join(age_v0, by = "eid") %>%
    mutate(
      followup_year = followup_time / 365.25,
      Overall_acceleration_average_valid =
        ifelse(
          calibration == 1 & good_wear_time == 1,
          as.numeric(Overall_acceleration_average),
          NA_real_
        )
    )
  
  # 循环生活方式
  for (life_var in lifestyle_vars) {
    message("Lifestyle: ", life_var)
    
    d[[life_var]] <- as.numeric(as.character(d[[life_var]]))
    
    dat <- d %>%
      dplyr::select(age_gap_adjusted, all_of(life_var), all_of(covars)) %>%
      filter(complete.cases(.))
    dat$life <- dat[[life_var]]
    
    # datadist
    dd <- datadist(dat)
    options(datadist = "dd")
    
    # k = 3-5
    fits <- list()
    aic <- data.frame()
    
    for (k in 3:5) {
      formula <- as.formula(paste0("age_gap_adjusted ~ rcs(life, ", k, ") + ", 
                                   paste(covars, collapse = " + ")))
      
      fit <- ols(formula, data = dat, x = TRUE, y = TRUE)
      fits[[paste0("k", k)]] <- fit
      
      aic <- rbind(
        aic,
        data.frame(
          outcome = outcome,
          lifestyle = life_var,
          k = k,
          N = nrow(dat),
          AIC = as.numeric(AIC(fit)[1])
        )
      )
    }
    
    # Best k
    best_k <- aic$k[which.min(aic$AIC)]
    aic$best_k <- best_k
    aic$delta_AIC <- aic$AIC - min(aic$AIC)
    
    AIC_result[[paste(outcome, life_var, sep = "__")]] <- aic
    
    # Prediction
    for (k in 3:5) {
      fit <- fits[[paste0("k", k)]]
      
      p <- as.data.frame(rms::Predict(fit, life, np = 200, conf.int = 0.95, ref.zero=TRUE))
      p <- p %>% mutate(k = k, outcome = outcome, lifestyle = life_var, best_k = best_k)
      
      pred_result[[paste(outcome, life_var, k, sep = "__")]] <- p
    }
  }
}


# 合并并保存 AIC

AIC_result <- bind_rows(AIC_result)
fwrite(AIC_result, file.path(out_dir, "RCS_AIC.csv"))

# 合并并保存 prediction
pred_result <- bind_rows(pred_result)
fwrite(pred_result, file.path(out_dir, "RCS_predictions.csv"))


# 绘图
pred_result <- fread(file.path(out_dir, "RCS_predictions.csv"))
head(pred_result)

out_dir <- "result/F2_生活方式RCS/plot"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

outcomes <- unique(pred_result$outcome)

for (life_var in unique(pred_result$lifestyle)) {
  d <- pred_result %>% filter(lifestyle == life_var)
  
  plots <- list()
  
  # 每个 outcome
  for (outcome_i in outcomes) {
    # k = 3,4,5
    for (k_i in 3:5) {
      pdat <- d %>% filter(outcome == outcome_i, k == k_i)
      
      p <- ggplot(pdat, aes(x = life, y = yhat)) +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(linewidth = 0.9) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        theme_bw() +
        theme(plot.title = element_text(size = 11, face = "bold")) +
        labs(title = paste0(outcome_i, "\nk = ", k_i), x = life_var, y = "V2 age acceleration")
      
      plots[[paste(outcome_i, paste0("k", k_i), sep = "__")]] <- p
    }
    
    # Best k
    best_k_i <- unique(d %>% filter(outcome == outcome_i) %>% pull(best_k))
    
    pdat <- d %>% filter(outcome == outcome_i, k == best_k_i)
    
    if (nrow(pdat) > 0) {
      p <- ggplot(pdat, aes(x = life, y = yhat)) +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
        geom_line(linewidth = 1) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        theme_bw() +
        theme(plot.title = element_text(size = 11, face = "bold")) +
        labs(title = paste0(outcome_i, "\nBest k = ", best_k_i), x = life_var, y = "Age acceleration")
      
      plots[[paste(outcome_i, "best", sep = "__")]] <- p
    }
  }
  # 按 outcome × k 排列成 3 × 4
  final_plot <- wrap_plots(plots, ncol = 4) +
    plot_annotation(title = life_var,
                    theme = theme(plot.title = element_text(size = 16, face = "bold")))
  
  # 保存
  ggsave(file.path(out_dir, paste0(life_var, "_RCS.png")), final_plot, width = 16, height = 11, dpi = 300)
}
