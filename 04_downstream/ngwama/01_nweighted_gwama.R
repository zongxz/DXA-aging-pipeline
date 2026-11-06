# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# N-weighted 多变量 GWAMA（用 LDSC 截距矩阵校正样本重叠）
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------

library(data.table)

# 依赖第三方函数 multivariate_GWAMA()（N-weighted GWAMA，v1.2.6）
# 源码: https://github.com/baselmans/multivariate_GWAMA
# 注意：必须用 Test_Data/N_weighted_GWAMA.function.1_2_6.R（函数名小写、含 h2 参数）；
#       同仓库 Software_used_in_paper/ 下的版本是另一个函数（Multivariate_GWAMA，无 h2），不可混用。
source("/path/to/tools/N_weighted_GWAMA.function.1_2_6.R")

knee <- fread("GWAS/gwas_result/single_right_knee.txt.gz")
hip  <- fread("GWAS/gwas_result/single_right_hip.txt.gz")

colnames(knee)
colnames(hip)

# Prepare NGWAMA input

knee_ngwama <- knee[, .(
  SNPID = SNP,
  CHR   = as.integer(CHR),
  BP    = as.integer(POS),
  EA    = A1,
  OA    = A2,
  EAF   = as.numeric(AF1),
  N     = as.numeric(N),
  Z     = as.numeric(BETA / SE),
  PVAL  = as.numeric(P)
)]

hip_ngwama <- hip[, .(
  SNPID = SNP,
  CHR   = as.integer(CHR),
  BP    = as.integer(POS),
  EA    = A1,
  OA    = A2,
  EAF   = as.numeric(AF1),
  N     = as.numeric(N),
  Z     = as.numeric(BETA / SE),
  PVAL  = as.numeric(P)
)]


# Convert to ordinary data.frame
knee_ngwama <- as.data.frame(knee_ngwama)
hip_ngwama  <- as.data.frame(hip_ngwama)

fwrite(knee_ngwama, "GWAS/gwas_result/single_right_knee_NGWAMA.txt.gz", sep = "\t")
fwrite(hip_ngwama, "GWAS/gwas_result/single_right_hip_NGWAMA.txt.gz", sep = "\t")


# Nweighted GWAMA ====
# 读取数据

dat<-vector("list",2)

dat[[1]]<-fread("GWAS/gwas_result/single_right_knee_NGWAMA.txt.gz",showProgress=F,data.table=F)
dat[[2]]<-fread("GWAS/gwas_result/single_right_hip_NGWAMA.txt.gz",showProgress=F,data.table=F)




# SNP遗传力（SNP heritability）
# 0.1017：右膝年龄加速（single_right_knee）的SNP遗传力
# 来源：LDSC中 single_right_knee 的 "Total Observed scale h2"
h2list <- c(
  0.1017,  # single_right_knee：LDSC Total Observed scale h2
  0.1247   # single_right_hip：LDSC Total Observed scale h2
)

# 跨性状LDSC截距矩阵（Cross-Trait Intercept, CTI）
# 对角线：各性状自身的LDSC截距（h2_int）
# 非对角线：两个性状之间的跨性状LDSC截距（gcov_int）
#
# 来源：
# 1.0005：single_right_knee的LDSC "h2_int"
# 1.0077：single_right_hip的LDSC "h2_int"
# 0.2286：single_right_knee × single_right_hip的LDSC "gcov_int"
#
# CTI矩阵的性状顺序必须与GWAS summary statistics列表 dat 中的顺序完全一致：
# dat[[1]] = single_right_knee
# dat[[2]] = single_right_hip
CTI <- matrix(
  c(
    1.0005, 0.2286,
    0.2286, 1.0077
  ),
  nrow = 2,
  byrow = TRUE
)

rownames(CTI) <- colnames(CTI) <- c(
  "single_right_knee",
  "single_right_hip"
)

# Check
CTI
h2list

multivariate_GWAMA(
  x = dat,
  cov_Z = CTI,
  h2 = h2list,
  out = "GWAS/NGWAMA/result",
  name = "knee_hip",
  output_gz = TRUE,
  check_columns = TRUE
)


result <- fread('GWAS/NGWAMA/result/knee_hip.N_weighted_GWAMA.results.txt.gz')
head(result)


nrow(knee_ngwama)
nrow(hip_ngwama)

length(unique(knee_ngwama$SNPID))
length(unique(hip_ngwama$SNPID))

length(intersect(
  knee_ngwama$SNPID,
  hip_ngwama$SNPID
))

length(union(
  knee_ngwama$SNPID,
  hip_ngwama$SNPID
))


# Model Averaging GWAMA ====
packages <- c(
  "sp",
  "raster",
  "RcppArmadillo",
  "unmarked",
  "VGAM",
  "AICcmodavg",
  "metafor"
)

sapply(packages, requireNamespace, quietly = TRUE)
install.packages(c('raster'))





# Result ====
N_weighted_result <- fread('GWAS/NGWAMA/result/knee_hip.N_weighted_GWAMA.results.txt.gz')
MA_result <- fread('GWAS/NGWAMA/result/knee_hip_MA_GWAMA.txt.gz')

head(N_weighted_result)
head(MA_result)




