# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# LAVA 结果整理：FDR 校正、局部 h2/遗传相关作图、显著区域导出
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
library(LAVA)
library(dplyr)
library(ggplot2)
library(foreach)
library(doParallel)
library(data.table)

# LAVA ====
# sample.overlap.txt
# 首先生成 all.rg
#
# FILES=($(ls *_rg.log))        # assuming the output format is [phenotype]_rg.log
# N=$(echo ${#FILES[@]})          # and that all combinations of penotypes have been analysed
#
#   for I in ${FILES[@]}; do
#   PHEN=$(echo $I | sed 's/_rg\.log//')
#
#   # subset log files to relevant output
#   tail -n$(($N+4)) $I | head -$((N+1)) > $PHEN.rg     # (adapt as necessary)
#
#   # add to single data set
#   if [[ $I == ${FILES[0]} ]]; then
#   cat $PHEN.rg > all.rg       # only including the header for the first phenotypes
#   else
#     cat $PHEN.rg | sed '1d' >> all.rg
#   fi
#   done
# 修改all.rg文件内容


scor = fread("/path/to/project/GWAS/LDSC/all.rg")              # read in
scor = scor[, c("p1", "p2", "gcov_int")]             # retain key headers
scor$p1 = gsub(".sumstats.gz", "", scor$p1)   # assuming the munged files have format [phenotype]_munge.sumstats.gz
scor$p2 = gsub(".sumstats.gz", "", scor$p2)   # (adapt as necessary)

phen <- unique(c(scor$p1, scor$p2))

mat <- matrix(
  NA_real_,
  nrow = length(phen),
  ncol = length(phen),
  dimnames = list(phen, phen)
)

# 先填入已有的 LDSC gcov_int
for (k in seq_len(nrow(scor))) {
  mat[scor$p1[k], scor$p2[k]] <- scor$gcov_int[k]
}
mat
# 自动用已有方向补齐另一方向
idx <- is.na(mat) & !is.na(t(mat))
mat[idx] <- t(mat)[idx]
mat

mat = round(cov2cor(mat), 5)                       # standardise
write.table(mat, "GWAS/LAVA/sample.overlap.txt", quote = F)

sample_cor_kh <- mat[2:3, 2:3]
sample_cor_kh
write.table(sample_cor_kh, file = "GWAS/LAVA/sample.overlap.knee.hip.txt", quote = F)

# LAVA 三个trait
input = process.input(
  input.info.file = "GWAS/LAVA/input.info.txt",
  sample.overlap.file = "GWAS/LAVA/sample.overlap.txt",
  ref.prefix = "/path/to/ref/g1000_eur/g1000_eur",
  phenos = c("double_right", "single_right_knee", "single_right_hip")
)       # subset of phenotypes listed in the input info file that we want to process


# LAVA 膝+髋
input = process.input(
  input.info.file = "GWAS/LAVA/input.info.knee.hip.txt",
  sample.overlap.file = "GWAS/LAVA/sample.overlap.knee.hip.txt",
  ref.prefix = "/path/to/ref/g1000_eur/g1000_eur",
  phenos = c("single_right_knee", "single_right_hip")
)



# LAVA Knee vs Knee − Hip
input <- process.input(
  input.info.file = "GWAS/LAVA/input.info.knee_double.txt",
  sample.overlap.file = "GWAS/LAVA/sample.overlap.knee_double.txt",
  ref.prefix = "/path/to/ref/g1000_eur/g1000_eur",
  phenos = c("single_right_knee", "double_right")
)

# LAVA Hip vs Knee − Hip
input <- process.input(
  input.info.file = "GWAS/LAVA/input.info.hip_double.txt",
  sample.overlap.file = "GWAS/LAVA/sample.overlap.hip_double.txt",
  ref.prefix = "/path/to/ref/g1000_eur/g1000_eur",
  phenos = c("single_right_hip", "double_right")
)


loci = read.loci("GWAS/LAVA/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")

locus <- process.locus(loci[1, ], input, drop.failed = FALSE)


registerDoParallel(cores = 20)

result = foreach(i = 1:nrow(loci), .packages = "LAVA") %dopar% {
  
  locus = process.locus(loci[i, ], input, drop.failed = FALSE)
  
  result = list()
  
  if (!is.null(locus)) {
    
    res_univ = cbind(loci = i, run.univ(locus))
    result[[1]] = res_univ
    
    if (sum(locus$failed) == 0) {
      
      res_bivar = tryCatch(
        cbind(loci = i, run.bivar(locus)),
        error = function(e) {
          message("Locus ", i, " bivariate failed: ",
                  conditionMessage(e))
          NULL
        }
      )
      
      if (!is.null(res_bivar)) {
        result[[2]] = res_bivar
      }
    }
  }
  
  result
}
# 注意：上面 4 次 process.input() 会依次覆盖 input，此处实际生效的是最后一个（hip vs double）。
# 因此保存名必须与之一致，否则会用 hip 的结果覆盖掉真正的 knee_double 结果。
# 每个性状对建议各跑一次 01_run_pair.R，只保留一次 process.input()，避免此陷阱。
save(result, file = 'GWAS/LAVA/lava.hip.double.Rdata')

# a====
# 结果====
# 1. 参数
base_dir <- "GWAS/LAVA"
phenotypes <- c("single_right_knee", "single_right_hip")
phenotypes <- c("single_right_knee", "single_right_hip", 'double_right')
phenotype_colors <- c(single_right_knee = "#006affff",
                      single_right_hip  = "#1dcc2fff",
                      single_right_hip  = "#FD3E5Cff")

# Knee VS Hip ====
# 2. 读取并整理 LAVA 结果
# 处理结果
library(LAVA)
loci <- read.loci("GWAS/LAVA/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")

load("GWAS/LAVA/lava.knee.hip.Rdata")

# univariate 
res_univ <- lapply(result, function(x) {
  if (length(x) > 0) { x[[1]] }
})

res_univ <- do.call(rbind, res_univ)

# 查看两个 phenotype 各有多少结果
table(res_univ$phen)

# 分别提取 knee 和 hip
res_univ_knee <- res_univ[res_univ$phen == "single_right_knee", ]
res_univ_hip <- res_univ[res_univ$phen == "single_right_hip", ]
# 按 loci 合并
res_univ <- merge(res_univ_knee[, -2], res_univ_hip[, -2], by = "loci", all = TRUE)

# 重命名
res_univ <- dplyr::rename(
  res_univ,
  p.phen.knee = p.x,
  h2.phen.knee = h2.obs.x,
  p.phen.hip = p.y,
  h2.phen.hip = h2.obs.y
)


# 加入 locus 的染色体和位置
res_univ <- merge(loci, res_univ, by.x = "LOC", by.y = "loci")

# bivariate 
res_bivar <- lapply(result, function(x) {
  if (length(x) > 1) { x[[2]] }
})

# 去掉没有 bivariate 结果的 locus
res_bivar <- res_bivar[!sapply(res_bivar, is.null)]

res_bivar <- do.call(rbind, res_bivar)

# 重命名
res_bivar <- dplyr::rename(res_bivar, p.rho = p)

# combine  
res_combine <- merge(res_univ, res_bivar[, c(1, 4:6, 10)], 
                     by.x = "LOC", by.y = "loci", all = TRUE)

# FDR correction for local genetic correlation
res_combine$p_knee_fdr <- p.adjust(res_combine$p.phen.knee, method = "BH")
res_combine$p_hip_fdr <- p.adjust(res_combine$p.phen.hip, method = "BH")
res_combine$p.rho_fdr <- p.adjust(res_combine$p.rho, method = "BH")

write.csv(res_combine,'GWAS/LAVA/result/res_combine_knee_hip.csv',row.names = F)
# region plot: Genetic Correlation
res_combine <- fread('GWAS/LAVA/result/res_combine_knee_hip.csv')
d1 <- res_combine %>%
  dplyr::filter(p_knee_fdr < 0.05 & p_hip_fdr < 0.05) %>%
  dplyr::select(LOC, CHR, START, STOP, rho, p.rho_fdr)


# loci with significant local genetic correlation
loc <- d1$LOC[d1$p.rho_fdr < 0.05]
loc <- loc[!is.na(loc)]

# labels
loc_labl <- d1[d1$LOC %in% loc, ]
loc_labl <- paste0(loc_labl$CHR, ":", loc_labl$START, "-", loc_labl$STOP)

d1 <- d1 %>% dplyr::select(-p.rho_fdr, -STOP)
d1
summary(d1$rho)
max(d1$rho, na.rm = TRUE)

nrow(d1)
sig_rg <- d1 %>%
  dplyr::filter(!is.na(rho), LOC %in% loc)

# 红色区域数量及相关方向
nrow(sig_rg)
table(
  direction = ifelse(
    sig_rg$rho > 0, "Positive",
    ifelse(sig_rg$rho < 0, "Negative", "Zero")
  )
)



CMplot(d1, plot.type = "m", pch = NA, type = "h", col = "grey", 
       cex = c(0.5, 0.5), signal.cex = c(0.5, 0.5), LOG10 = F, 
       highlight = loc, highlight.type = "h", highlight.col = "red", 
       highlight.pch = NA, highlight.cex = 0.8, 
       #highlight.text = loc_labl,  highlight.text.cex = 0.7,highlight.text.font = 2,
       file = "pdf",
       file.name = "hip_knee_cor",  ylab = "Genetic Correlation", height = 5, width = 17 )


# Knee VS Knee - Hip ====
# region plot: Genetic Correlation
# 2. 读取并整理 LAVA 结果
# 处理结果
library(LAVA)
loci <- read.loci("GWAS/LAVA/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")

load("GWAS/LAVA/lava.knee.double.Rdata")

# univariate 
res_univ <- lapply(result, function(x) {
  if (length(x) > 0) { x[[1]] }
})

res_univ <- do.call(rbind, res_univ)

# 查看两个 phenotype 各有多少结果
table(res_univ$phen)

# 分别提取 knee 和 double
res_univ_knee <- res_univ[res_univ$phen == "single_right_knee", ]
res_univ_double <- res_univ[res_univ$phen == "double_right", ]
# 按 loci 合并
res_univ <- merge(res_univ_knee[, -2], res_univ_double[, -2], by = "loci", all = TRUE)

# 重命名
res_univ <- dplyr::rename(
  res_univ,
  p.phen.knee = p.x,
  h2.phen.knee = h2.obs.x,
  p.phen.double = p.y,
  h2.phen.double = h2.obs.y
)


# 加入 locus 的染色体和位置
res_univ <- merge(loci, res_univ, by.x = "LOC", by.y = "loci")

# bivariate 
res_bivar <- lapply(result, function(x) {
  if (length(x) > 1) { x[[2]] }
})

# 去掉没有 bivariate 结果的 locus
res_bivar <- res_bivar[!sapply(res_bivar, is.null)]

res_bivar <- do.call(rbind, res_bivar)

# 重命名
res_bivar <- dplyr::rename(res_bivar, p.rho = p)

# combine  
res_combine <- merge(res_univ, res_bivar[, c(1, 4:6, 10)], 
                     by.x = "LOC", by.y = "loci", all = TRUE)

# FDR correction for local genetic correlation
res_combine$p_knee_fdr <- p.adjust(res_combine$p.phen.knee, method = "BH")
res_combine$p_double_fdr <- p.adjust(res_combine$p.phen.double, method = "BH")
res_combine$p.rho_fdr <- p.adjust(res_combine$p.rho, method = "BH")

write.csv(res_combine,'GWAS/LAVA/result/res_combine_knee_double.csv',row.names = F)
# region plot: Genetic Correlation
res_combine <- fread('GWAS/LAVA/result/res_combine_knee_double.csv')
d1 <- res_combine %>%
  dplyr::filter(p_knee_fdr < 0.05 & p_double_fdr < 0.05) %>%
  dplyr::select(LOC, CHR, START, STOP, rho, p.rho_fdr)


# loci with significant local genetic correlation
loc <- d1$LOC[d1$p.rho_fdr < 0.05]
loc <- loc[!is.na(loc)]

# labels
loc_labl <- d1[d1$LOC %in% loc, ]
loc_labl <- paste0(loc_labl$CHR, ":", loc_labl$START, "-", loc_labl$STOP)

d1 <- d1 %>% dplyr::select(-p.rho_fdr, -STOP)
d1
summary(d1$rho)
max(d1$rho, na.rm = TRUE)


nrow(d1)
sig_rg <- d1 %>%
  dplyr::filter(!is.na(rho), LOC %in% loc)

# 红色区域数量及相关方向
nrow(sig_rg)
table(
  direction = ifelse(
    sig_rg$rho > 0, "Positive",
    ifelse(sig_rg$rho < 0, "Negative", "Zero")
  )
)



CMplot(d1, plot.type = "m", pch = NA, type = "h", col = "grey", 
       cex = c(0.5, 0.5), signal.cex = c(0.5, 0.5), LOG10 = F, 
       highlight = loc, highlight.type = "h", highlight.col = "red", 
       highlight.pch = NA, highlight.cex = 0.8, 
       #highlight.text = loc_labl,  highlight.text.cex = 0.7,highlight.text.font = 2,
       file = "pdf",ylim = c(0, 1),
       file.name = "double_knee_cor",  ylab = "Genetic Correlation", height = 5, width = 17 )


# Hip VS Knee - Hip ====
# 2. 读取并整理 LAVA 结果
# 处理结果 
library(LAVA)
loci <- read.loci("GWAS/LAVA/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")

load("GWAS/LAVA/lava.hip.double.Rdata")

# univariate 
res_univ <- lapply(result, function(x) {
  if (length(x) > 0) { x[[1]] }
})

res_univ <- do.call(rbind, res_univ)

# 查看两个 phenotype 各有多少结果
table(res_univ$phen)

# 分别提取 hip 和 double
res_univ_hip <- res_univ[res_univ$phen == "single_right_hip", ]
res_univ_double <- res_univ[res_univ$phen == "double_right", ]
# 按 loci 合并
res_univ <- merge(res_univ_hip[, -2], res_univ_double[, -2], by = "loci", all = TRUE)

# 重命名
res_univ <- dplyr::rename(
  res_univ,
  p.phen.hip = p.x,
  h2.phen.hip = h2.obs.x,
  p.phen.double = p.y,
  h2.phen.double = h2.obs.y
)


# 加入 locus 的染色体和位置
res_univ <- merge(loci, res_univ, by.x = "LOC", by.y = "loci")

# bivariate 
res_bivar <- lapply(result, function(x) {
  if (length(x) > 1) { x[[2]] }
})

# 去掉没有 bivariate 结果的 locus
res_bivar <- res_bivar[!sapply(res_bivar, is.null)]

res_bivar <- do.call(rbind, res_bivar)

# 重命名
res_bivar <- dplyr::rename(res_bivar, p.rho = p)

# combine  
res_combine <- merge(res_univ, res_bivar[, c(1, 4:6, 10)], 
                     by.x = "LOC", by.y = "loci", all = TRUE)

# FDR correction for local genetic correlation
res_combine$p_hip_fdr <- p.adjust(res_combine$p.phen.hip, method = "BH")
res_combine$p_double_fdr <- p.adjust(res_combine$p.phen.double, method = "BH")
res_combine$p.rho_fdr <- p.adjust(res_combine$p.rho, method = "BH")

write.csv(res_combine,'GWAS/LAVA/result/res_combine_hip_double.csv',row.names = F)
# region plot: Genetic Correlation
res_combine <- fread('GWAS/LAVA/result/res_combine_hip_double.csv')
d1 <- res_combine %>%
  dplyr::filter(p_hip_fdr < 0.05 & p_double_fdr < 0.05) %>%
  dplyr::select(LOC, CHR, START, STOP, rho, p.rho_fdr)


# loci with significant local genetic correlation
loc <- d1$LOC[d1$p.rho_fdr < 0.05]
loc <- loc[!is.na(loc)]

# labels
loc_labl <- d1[d1$LOC %in% loc, ]
loc_labl <- paste0(loc_labl$CHR, ":", loc_labl$START, "-", loc_labl$STOP)

d1 <- d1 %>% dplyr::select(-p.rho_fdr, -STOP)
d1
summary(d1$rho)
max(d1$rho, na.rm = TRUE)


nrow(d1)
sig_rg <- d1 %>%
  dplyr::filter(!is.na(rho), LOC %in% loc)

# 红色区域数量及相关方向
nrow(sig_rg)
table(
  direction = ifelse(
    sig_rg$rho > 0, "Positive",
    ifelse(sig_rg$rho < 0, "Negative", "Zero")
  )
)



CMplot(d1, plot.type = "m", pch = NA, type = "h", col = "grey", 
       cex = c(0.5, 0.5), signal.cex = c(0.5, 0.5), LOG10 = F, 
       highlight = loc, highlight.type = "h", highlight.col = "red", 
       highlight.pch = NA, highlight.cex = 0.8, 
       #highlight.text = loc_labl,  highlight.text.cex = 0.7,highlight.text.font = 2,
       file = "pdf",
       file.name = "double_hip_cor",  ylab = "Genetic Correlation", height = 5, width = 17 )

# local heritability ====
# Hip local heritability 
d_hip <- res_combine %>% dplyr::select(LOC, CHR, START, h2.phen.hip)

# Hip FDR < 0.05 的 loci
loc_hip <- res_combine$LOC[res_combine$p_hip_fdr < 0.05]

loc_hip <- loc_hip[!is.na(loc_hip)]

# Hip 显著 locus 标签
loc_labl_hip <- res_combine[res_combine$LOC %in% loc_hip, ]

loc_labl_hip <- paste0(loc_labl_hip$CHR, ":", loc_labl_hip$START, "-", loc_labl_hip$STOP)

CMplot(
  d_hip,
  plot.type = "m",
  type = "h",
  pch = NA,
  col = "grey",
  threshold.col = c("grey", "black"),
  threshold.lty = c(1, 1),
  threshold.lwd = c(1, 1),
  amplify = TRUE,
  cex = c(0.5, 0.5),
  signal.cex = c(0.5, 0.5),
  LOG10 = FALSE,
  highlight = loc_hip,
  highlight.type = "h",
  highlight.col = "#1dcc2fff",
  highlight.pch = NA,
  #highlight.text = loc_labl_hip,
  highlight.cex = 1.2,
  file = "pdf",
  file.name = "hip",
  ylab = "Heritability",
  height = 5,
  width = 17
)


# Knee local heritability 

d_knee <- res_combine %>%
  dplyr::select(LOC, CHR, START, h2.phen.knee)

# Knee FDR < 0.05 的 loci
loc_knee <- res_combine$LOC[res_combine$p_knee_fdr < 0.05]

loc_knee <- loc_knee[!is.na(loc_knee)]

# Knee 显著 locus 标签
loc_labl_knee <- res_combine[res_combine$LOC %in% loc_knee, ]

loc_labl_knee <- paste0(loc_labl_knee$CHR,
                        ":",
                        loc_labl_knee$START,
                        "-",
                        loc_labl_knee$STOP)

CMplot(
  d_knee,
  plot.type = "m",
  type = "h",
  pch = NA,
  col = "grey",
  
  threshold.col = c("grey", "black"),
  threshold.lty = c(1, 1),
  threshold.lwd = c(1, 1),
  amplify = TRUE,
  
  cex = c(0.5, 0.5),
  signal.cex = c(0.5, 0.5),
  LOG10 = FALSE,
  
  highlight = loc_knee,
  highlight.type = "h",
  highlight.col = "#006affff",
  highlight.pch = NA,
  #highlight.text = loc_labl_knee,
  highlight.cex = 1.2,
  
  file = "pdf",
  file.name = "knee",
  ylab = "Heritability",
  height = 5,
  width = 17
)

# Knee-Hip local heritability 
d_double <- res_combine %>% dplyr::select(LOC, CHR, START, h2.phen.double)

# Knee-Hip FDR < 0.05 的 loci
loc_double <- res_combine$LOC[res_combine$p_double_fdr < 0.05]

loc_double <- loc_double[!is.na(loc_double)]
# Knee-Hip 显著 locus 标签
loc_labl_double <- res_combine[res_combine$LOC %in% loc_double, ]

loc_labl_double <- paste0(loc_labl_double$CHR, ":", loc_labl_double$START, "-", loc_labl_double$STOP)

CMplot(
  d_double,
  plot.type = "m",
  type = "h",
  pch = NA,
  col = "grey",
  threshold.col = c("grey", "black"),
  threshold.lty = c(1, 1),
  threshold.lwd = c(1, 1),
  amplify = TRUE,
  cex = c(0.5, 0.5),
  signal.cex = c(0.5, 0.5),
  LOG10 = FALSE,
  highlight = loc_double,
  highlight.type = "h",
  highlight.col = "#FD3E5Cff",
  highlight.pch = NA,
  #highlight.text = loc_labl_hip,
  highlight.cex = 1.2,
  file = "pdf",
  file.name = "double",
  ylab = "Heritability",
  height = 5,
  width = 17
)


# 提取数据====
## Knee VS Hip ====
res_combine <- read.csv("GWAS/LAVA/result/res_combine_knee_hip.csv")

# 提取 34 个显著局部遗传相关区域
lava_sig <- res_combine %>%
  filter(p_knee_fdr < 0.05, p_hip_fdr < 0.05, p.rho_fdr < 0.05, !is.na(rho)) %>%
  select(
    CHR, START, STOP, h2.phen.knee, p.phen.knee, h2.phen.hip, 
    p.phen.hip, rho, p.rho, p_knee_fdr, p_hip_fdr, p.rho_fdr
  ) %>%
  rename(
    Chromosome = CHR,
    Start_position = START,
    End_position = STOP,
    Knee_Heritability = h2.phen.knee,
    Knee_Pvalue = p.phen.knee,
    Hip_Heritability = h2.phen.hip,
    Hip_Pvalue = p.phen.hip,
    Genetic_Correlation = rho,
    Genetic_Correlation_Pvalue = p.rho
  ) %>%
  arrange(Chromosome, Start_position)

# 检查提取数量
nrow(lava_sig)

# 查看结果
print(lava_sig)

# 保存结果
write.csv(lava_sig,"GWAS/LAVA/result/LAVA_significant_knee_hip_34_regions.csv",row.names = F)

## Knee VS Knee - double ====
res_combine <- read.csv("GWAS/LAVA/result/res_combine_knee_double.csv")

# 提取 显著局部遗传相关区域
lava_sig <- res_combine %>%
  filter(p_knee_fdr < 0.05, p_double_fdr < 0.05, p.rho_fdr < 0.05, !is.na(rho)) %>%
  select(
    CHR, START, STOP, h2.phen.knee, p.phen.knee, h2.phen.double, 
    p.phen.double, rho, p.rho, p_knee_fdr, p_double_fdr, p.rho_fdr
  ) %>%
  rename(
    Chromosome = CHR,
    Start_position = START,
    End_position = STOP,
    Knee_Heritability = h2.phen.knee,
    Knee_Pvalue = p.phen.knee,
    double_Heritability = h2.phen.double,
    double_Pvalue = p.phen.double,
    Genetic_Correlation = rho,
    Genetic_Correlation_Pvalue = p.rho
  ) %>%
  arrange(Chromosome, Start_position)

# 检查提取数量
nrow(lava_sig)

# 查看结果
print(lava_sig)

# 保存结果
write.csv(lava_sig,"GWAS/LAVA/result/LAVA_significant_knee_double_207_regions.csv",row.names = F)
## Hip VS hip - double ====
res_combine <- read.csv("GWAS/LAVA/result/res_combine_hip_double.csv")

# 提取 显著局部遗传相关区域
lava_sig <- res_combine %>%
  filter(p_hip_fdr < 0.05, p_double_fdr < 0.05, p.rho_fdr < 0.05, !is.na(rho)) %>%
  select(
    CHR, START, STOP, h2.phen.hip, p.phen.hip, h2.phen.double, 
    p.phen.double, rho, p.rho, p_hip_fdr, p_double_fdr, p.rho_fdr
  ) %>%
  rename(
    Chromosome = CHR,
    Start_position = START,
    End_position = STOP,
    hip_Heritability = h2.phen.hip,
    hip_Pvalue = p.phen.hip,
    double_Heritability = h2.phen.double,
    double_Pvalue = p.phen.double,
    Genetic_Correlation = rho,
    Genetic_Correlation_Pvalue = p.rho
  ) %>%
  arrange(Chromosome, Start_position)

# 检查提取数量
nrow(lava_sig)

# 查看结果
print(lava_sig)

# 保存结果
write.csv(lava_sig,"GWAS/LAVA/result/LAVA_significant_hip_double_301_regions.csv",row.names = F)
