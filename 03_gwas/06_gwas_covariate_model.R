# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 03_gwas
# GWAS 前置：表型/协变量整理、样本交集核对与模型设定
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
library(data.table)
# GWAS数据处理=====

eid <- fread('eye/data/eid_match_ukb.csv')
colnames(eid) = c('eid1','eid')


# 批次数据
batch <- fread('/path/to/ref/ukb/Genotyping process and sample QC.csv')
batch <- batch[,1:2]
colnames(batch) <- c('eid', 'batch')
batch[batch == ""] <- NA

# PC数据
pc <- fread('/path/to/ref/ukb/PC40.txt')
colnames(pc)[colnames(pc) == 'f.eid'] <- 'eid' 
# 基因型性别
g_sex <- fread('/path/to/ref/ukb/genetic_sex.csv')
colnames(g_sex) <- c('eid', 'sex')

# covar (基因型性别和批次)
covar <- fread('/path/to/ref/ukb/covar.txt')
colnames(covar)[colnames(covar) == 'V1'] <- 'eid1' 
colnames(covar)[colnames(covar) == 'V2'] <- 'eid' 

# qcovar (PC和年龄)
qcovar <- fread('/path/to/ref/ukb/qcovar.txt')
colnames(qcovar)[colnames(covar) == 'V1'] <- 'eid1' 
colnames(qcovar)[colnames(covar) == 'V2'] <- 'eid' 


# 提取协变量
data.cov <- inner_join(g_sex, batch, by='eid')
data.cov$batch <- sub("_.*", "", data.cov$batch)

data.cov <- inner_join(eid, data.cov, by='eid')
data.cov$eid <- data.cov$eid1

data.cov <- inner_join(data.cov,pc[,c(1,3:22)])
# 使用 mutate 和 case_when 来处理
data.cov <- data.cov %>%
  mutate(sex = case_when(
    sex == "Female" ~ '0',
    sex == "Male" ~ '1',
    sex == "" ~ NA_character_,
    TRUE ~ sex  # 保持其他值不变
  )) %>%
  mutate(batch = case_when(
    batch == "Batch" ~ '0',
    batch == "UKBiLEVEAX" ~ '1',
    TRUE ~ batch  # 保持其他值不变
  ))
write.csv(data.cov,'shared_project/data/gwas_cov.csv', row.names = F)





data.cov <- fread('shared_project/data/gwas_cov.csv')

# 年龄加速数据
double_right <- fread('data/gwas/double_right_oof.csv')
single_right_knee <- fread('data/gwas/single_right_knee_oof.csv')
single_right_hip <- fread('data/gwas/single_right_hip_oof.csv')

pred_age <- double_right

pred_age <- inner_join(eid,pred_age)
pred_age <- pred_age[,-2]
colnames(pred_age)[1] <- 'eid'

names(pred_age)

pred_age_gwas <- inner_join(pred_age[,c(1,8,9)], data.cov)
names(pred_age_gwas)

gwas_id <- fread('data/gwas/gwas_id.txt')
colnames(gwas_id)[1] <- 'eid'

pred_age_gwas <- pred_age_gwas %>% filter(eid %in% gwas_id$eid)

write.csv(pred_age_gwas, 'data/gwas/gwas_data_double_right.csv', row.names = F)
colnames(pred_age_gwas)

a = pred_age_gwas[, c(1,1,2)] # acc
b = pred_age_gwas[, c(1,1,5,6)] # 二分类协变量
c = pred_age_gwas[, c(1,1,3,7:26)] # 连续性协变量

colnames(c)
write.table(a, file = "/path/to/project/data/gwas/double_right.txt", sep = " ", 
            row.names = F, col.names = F, quote = TRUE)
write.table(b, file = "/path/to/project/data/gwas/covar.txt", sep = " ", 
            row.names = F, col.names = F, quote = TRUE)
write.table(c, file = "/path/to/project/data/gwas/qcovar.txt", sep = " ", 
            row.names = F, col.names = F, quote = TRUE)

# GWAS结果处理====
## double_right ====
double_right <- fread('GWAS/gwas_result/double_right.fastGWA')
double_right$MAF <- pmin(double_right$AF1, 1-double_right$AF1)
# fuma
write.table(double_right[,1:10], file="GWAS/gwas_result/double_right.txt",row.names = F,quote = F)

# scpagwas
head(double_right)
gwas_scpagwas <- double_right[,c(1,3,2,9,8,11)]
colnames(gwas_scpagwas)
colnames(gwas_scpagwas) <- c('chrom', 'pos', 'rsid', 'se', 'beta', 'maf')
write.table(gwas_scpagwas, file="GWAS/gwas_result/double_right_scpagwas.txt",row.names = F,quote = F)

# smr
results_auto = double_right %>% dplyr:: select(SNP,A1,A2,
                                       freq=AF1, b=BETA,
                                       se=SE, p=P, n=N)

fwrite(
  results_auto,
  "/path/to/project/GWAS/gwas_result/double_right_smr.txt.gz",
  sep = "\t"
)

# LD
ld_snp <- double_right[double_right$P<5e-8, c(2,10)]
colnames(ld_snp) <- c('rsid', 'pval')
library(ieugwasr)
bfile = "/path/to/ref/1kg.v3/EUR"  # path of reference file for LD computation
options(timeout = 600000000) 
clump_iv = ld_clump(dat = ld_snp, 
                    clump_kb = 500,
                    clump_r2 = 0.1,
                    clump_p = 1,
                    # opengwas_jwt=jwt,
                    plink_bin = '/path/to/tools/plink_19/plink',
                    bfile = bfile,
                    pop = "EUR")

# LAVA
fwrite(
  double_right[, .(SNP, A1, A2, N, BETA, P)],
  "/path/to/project/GWAS/gwas_result/double_right.lava.txt.gz",
  sep = "\t"
)

# PRS cs

fwrite(
  double_right[, .(SNP, A1, A2, BETA, P)],
  "/path/to/project/GWAS/gwas_result/double_right.prscs.txt",
  sep = "\t"
)

## single_right_knee ====
single_right_knee <- fread('GWAS/gwas_result/single_right_knee.fastGWA')
single_right_knee$MAF <- pmin(single_right_knee$AF1, 1-single_right_knee$AF1)

write.table(single_right_knee[,1:10], file="GWAS/gwas_result/single_right_knee.txt",row.names = F,quote = F)

# scpagwas
head(single_right_knee)
gwas_scpagwas <- single_right_knee[,c(1,3,2,9,8,11)]
colnames(gwas_scpagwas)
colnames(gwas_scpagwas) <- c('chrom', 'pos', 'rsid', 'se', 'beta', 'maf')
write.table(gwas_scpagwas, file="GWAS/gwas_result/single_right_knee_scpagwas.txt",row.names = F,quote = F)

# smr
results_auto = single_right_knee %>% dplyr:: select(SNP,A1,A2,
                                               freq=AF1, b=BETA,
                                               se=SE, p=P, n=N)
fwrite(
  results_auto,
  "/path/to/project/GWAS/gwas_result/single_right_knee_smr.txt.gz",
  sep = "\t"
)

# LD
ld_snp <- single_right_knee[single_right_knee$P<5e-8, c(2,10)]
colnames(ld_snp) <- c('rsid', 'pval')
library(ieugwasr)
bfile = "/path/to/ref/1kg.v3/EUR"  # path of reference file for LD computation
options(timeout = 600000000) 
clump_iv = ld_clump(dat = ld_snp, 
                    clump_kb = 500,
                    clump_r2 = 0.1,
                    clump_p = 1,
                    # opengwas_jwt=jwt,
                    plink_bin = '/path/to/tools/plink_19/plink',
                    bfile = bfile,
                    pop = "EUR"
)
# LAVA
fwrite(
  single_right_knee[, .(SNP, A1, A2, N, BETA, P)],
  "/path/to/project/GWAS/gwas_result/single_right_knee.lava.txt.gz",
  sep = "\t"
)
# PRS cs
fwrite(
  single_right_knee[, .(SNP, A1, A2, BETA, P)],
  "/path/to/project/GWAS/gwas_result/single_right_knee.prscs.txt",
  sep = "\t"
)

## single_right_hip ==== 
single_right_hip <- fread('GWAS/gwas_result/single_right_hip.fastGWA')
single_right_hip$MAF <- pmin(single_right_hip$AF1, 1-single_right_hip$AF1)

write.table(single_right_hip[,1:10], file="GWAS/gwas_result/single_right_hip.txt",row.names = F,quote = F)

# scpagwas
head(single_right_hip)
gwas_scpagwas <- single_right_hip[,c(1,3,2,9,8,11)]
colnames(gwas_scpagwas)
colnames(gwas_scpagwas) <- c('chrom', 'pos', 'rsid', 'se', 'beta', 'maf')
write.table(gwas_scpagwas, file="GWAS/gwas_result/single_right_hip_scpagwas.txt",row.names = F,quote = F)

# smr
results_auto = single_right_hip %>% dplyr:: select(SNP,A1,A2,
                                                    freq=AF1, b=BETA,
                                                    se=SE, p=P, n=N)
fwrite(
  results_auto,
  "/path/to/project/GWAS/gwas_result/single_right_hip_smr.txt.gz",
  sep = "\t"
)


# LD
ld_snp <- single_right_hip[single_right_hip$P<5e-8, c(2,10)]
colnames(ld_snp) <- c('rsid', 'pval')
library(ieugwasr)
bfile = "/path/to/ref/1kg.v3/EUR"  # path of reference file for LD computation
options(timeout = 600000000) 
clump_iv = ld_clump(dat = ld_snp, 
                    clump_kb = 500,
                    clump_r2 = 0.1,
                    clump_p = 1,
                    # opengwas_jwt=jwt,
                    plink_bin = '/path/to/tools/plink_19/plink',
                    bfile = bfile,
                    pop = "EUR"
)

# LAVA
fwrite(
  single_right_hip[, .(SNP, A1, A2, N, BETA, P)],
  "/path/to/project/GWAS/gwas_result/single_right_hip.lava.txt.gz",
  sep = "\t"
)
# PRS cs
fwrite(
  single_right_hip[, .(SNP, A1, A2, BETA, P)],
  "/path/to/project/GWAS/gwas_result/single_right_hip.prscs.txt",
  sep = "\t"
)


# CPASSOC ====
# LD
CPASSOC_result <- fread('GWAS/CPASSOC/CPASSOC_result.csv')
names(CPASSOC_result)
# P_Hom
ld_snp <- CPASSOC_result[CPASSOC_result$P_Hom<5e-8, c(1,5)]
colnames(ld_snp) <- c('rsid', 'pval')
length(unique(ld_snp$rsid))
length(ld_snp$rsid)
library(ieugwasr)
bfile = "/path/to/ref/1kg.v3/EUR"  # path of reference file for LD computation
options(timeout = 600000000) 
clump_iv = ld_clump(dat = ld_snp, 
                    clump_kb = 250,
                    clump_r2 = 0.6,
                    clump_p = 1,
                    # opengwas_jwt=jwt,
                    plink_bin = '/path/to/tools/plink_19/plink',
                    bfile = bfile,
                    pop = "EUR")


# fuma leadSNP
lead_snp <- fread('GWAS/CPASSOC/leadSNPs.txt')
ld_snp <- lead_snp$rsID

single_right_hip <- fread('GWAS/gwas_result/single_right_hip.fastGWA')
single_right_knee <- fread('GWAS/gwas_result/single_right_knee.fastGWA')
double_right <- fread('GWAS/gwas_result/double_right.fastGWA')
N_weighted_result <- fread('GWAS/NGWAMA/result/knee_hip.N_weighted_GWAMA.results.txt.gz')

hip <- single_right_hip[single_right_hip$SNP %in% ld_snp, ]
knee <-  single_right_knee[single_right_knee$SNP %in% ld_snp, ]
Knee_hip <-  double_right[double_right$SNP %in% ld_snp, ]
N_NGWAMA <-  N_weighted_result[N_weighted_result$SNPID %in% ld_snp, ]
CPASSOC <- lead_snp[,c(2,4,7)]



names(hip)[2] <- 'rsid'
names(knee)[2] <- 'rsid'
names(Knee_hip)[2] <- 'rsid'
names(N_NGWAMA)[2] <- 'rsid'
names(CPASSOC)[2] <- 'rsid'
names(CPASSOC)[3] <- 'P_Hom'
names(CPASSOC)


library(dplyr)
library(purrr)

# 给除 rsid 外的所有列加前缀
add_prefix <- function(df, prefix) {
  df %>% rename_with(~ paste0(prefix, "_", .x), .cols = -rsid)
}

# 分别加前缀
hip_pref    <- add_prefix(hip, "hip")
knee_pref   <- add_prefix(knee, "knee")
knee_hip_pref   <- add_prefix(Knee_hip, "knee_hip")
ngwama_pref <- add_prefix(N_NGWAMA, "N_NGWAMA")

# 按 rsid 合并所有数据框（全连接）
merged <- reduce(
  list(CPASSOC, hip_pref, knee_pref, knee_hip_pref,ngwama_pref),
  full_join,
  by = "rsid"
)

# 查看结果
head(merged)
write.csv(merged, 'GWAS/CPASSOC/lead_snp.csv',row.names = F)



CPASSOC_result <- fread('GWAS/CPASSOC/CPASSOC_result.csv')
names(CPASSOC_result)
lead_snp <- fread('GWAS/CPASSOC/leadSNPs.txt')
names(lead_snp) <- c("No", "GenomicLocus", "uniqID", "SNP", "CHR","POS", "P","nIndSigSNPs","IndSigSNPs")
gwas <- CPASSOC_result[CPASSOC_result$SNPID %in% lead_snp$SNP,]
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





