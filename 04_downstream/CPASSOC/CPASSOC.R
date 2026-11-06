data1 <- fread('/path/to/project/CPASSOC/CPASSOC/GWAS_data/single_right_hip.txt')
data2 <- fread('/path/to/project/CPASSOC/CPASSOC/GWAS_data/single_right_knee.txt')
#data3 <- fread('/home/zyw/total_data/sleeppoint/GWAS/results/19-21/fastGWA.fastGWA')

#构建SNP新的ID，避免重复rsID的情况
data1$variant_id <- paste(
  data1$CHR, data1$SNP, data1$POS, data1$A1, data1$A2,
  sep = ":"
)

data2$variant_id <- paste(
  data2$CHR, data2$SNP, data2$POS, data2$A1, data2$A2,
  sep = ":"
)

#data3$variant_id <- paste(
#  data3$CHR, data3$SNP, data3$POS, data3$A1, data3$A2,
#  sep = ":"
#)

head(data1[, c("CHR","SNP","POS","A1","A2","variant_id")])

#检查还有没有重复ID，为0才行
sum(duplicated(data1$variant_id))
sum(duplicated(data2$variant_id))
#sum(duplicated(data3$variant_id))

common_snp <- Reduce(
  intersect,
  list(
    data1$variant_id,
    data2$variant_id
    #data3$variant_id
  )
)

length(common_snp)

#只保留三者都有的SNP
d1 <- data1[data1$variant_id %in% common_snp, ]
d2 <- data2[data2$variant_id %in% common_snp, ]
#d3 <- data3[data3$variant_id %in% common_snp, ]

#结果为TRUE才行
setequal(d1$variant_id, d2$variant_id)
#setequal(d1$variant_id, d3$variant_id)

d1$Z1 <- d1$BETA / d1$SE
d2$Z2 <- d2$BETA / d2$SE
#d3$Z3 <- d3$BETA / d3$SE

d1 <- d1[, c("variant_id", "Z1")]
d2 <- d2[, c("variant_id", "Z2")]
#d3 <- d3[, c("variant_id", "Z3")]

Sumstat <- d1 %>%
  inner_join(d2, by = "variant_id")

#Sumstat <- d1 %>%
#  inner_join(d2, by = "variant_id") %>%
#  inner_join(d3, by = "variant_id")




#bash，将原始基因型数据rsID列重命名(这一段跑一次就行)
for chr in {1..22}; do

awk 'BEGIN{OFS="\t"} {
        print $1, $1":"$2":"$4":"$5":"$6, $3, $4, $5, $6
    }' \
/home/share/info_pass/${chr}.bim \
> /home/zyw/total_data/genotype/${chr}.bim

ln -sf /home/share/info_pass/${chr}.bed \
/home/zyw/total_data/genotype/${chr}.bed

ln -sf /home/share/info_pass/${chr}.fam \
/home/zyw/total_data/genotype/${chr}.fam

done


plink2 \
--bfile /home/zyw/total_data/genotype/1 \
--freq \
--out /home/zyw/total_data/genotype/test_chr1


library(data.table)

common_df <- data1[data1$variant_id %in% common_snp,
                   c("CHR", "variant_id")]

common_df <- unique(common_df)

for (chr in 1:22) {
  
  snps_chr <- common_df[
    CHR == chr,
    variant_id
  ]
  
  fwrite(
    data.table(variant_id = snps_chr),
    paste0(
      "/path/to/project/CPASSOC/CPASSOC/GWAS_data/common_chr",
      chr,
      ".txt"
    ),
    col.names = FALSE
  )
}


#bash
for chr in {1..22}; do

plink2 \
--bfile /home/zyw/total_data/genotype/${chr} \
--extract /path/to/project/CPASSOC/CPASSOC/GWAS_data/common_chr${chr}.txt \
--indep-pairwise 500kb 1 0.2 \
--out /path/to/project/CPASSOC/CPASSOC/LD_pruning/ld_chr${chr}

done

#合并22条染色体的prune.in
cat /path/to/project/CPASSOC/CPASSOC/LD_pruning/ld_chr{1..22}.prune.in \
> /path/to/project/CPASSOC/CPASSOC/LD_pruning/LD_independent_all.txt

head /path/to/project/CPASSOC/CPASSOC/LD_pruning/LD_independent_all.txt


#R
ld_snps <- fread(
  "/path/to/project/CPASSOC/CPASSOC/LD_pruning/LD_independent_all.txt",
  header = FALSE
)

colnames(ld_snps) <- "variant_id"

cor_data <- Sumstat %>%
  inner_join(ld_snps, by = "variant_id")

#去掉有缺失Z值的SNP
cor_data <- cor_data %>%
  filter(
    is.finite(Z1),
    is.finite(Z2)
    #is.finite(Z3)
  )

#排除Z大于1.96的
cor_data_null <- cor_data %>%
  filter(
    abs(Z1) <= 1.96,
    abs(Z2) <= 1.96
    #abs(Z3) <= 1.96
  )

#计算性状间的correlation矩阵
cor_matrix <- cor(
  cor_data_null[, c("Z1", "Z2")],
  use = "complete.obs",
  method = "pearson"
)

cor_matrix

write.table(
  cor_matrix,
  "/path/to/project/CPASSOC/CPASSOC/LD_pruning/CPASSOC_cor_matrix.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)



# CPASSOC analysis
# 1. Load CPASSOC functions
source("/home/zyw/total_data/CPASSOC/FunctionSet.R")

# 2. Prepare GWAS Z-score matrix
X <- as.matrix(
  Sumstat[, c("Z1", "Z2")]
)

# Make sure X is numeric
storage.mode(X) <- "numeric"


# 3. Sample sizes of the three GWAS
Samplesize <- c(
  48850,
  48850
)

# 4. Read correlation matrix
CorrMatrix <- as.matrix(
  read.table(
    "/path/to/project/CPASSOC/CPASSOC/LD_pruning/CPASSOC_cor_matrix.txt",
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
)

storage.mode(CorrMatrix) <- "numeric"


# 5. Check input data

cat("Number of SNPs:", nrow(X), "\n")
cat("Number of traits:", ncol(X), "\n")

cat("\nSample sizes:\n")
print(Samplesize)

cat("\nCorrelation matrix:\n")
print(CorrMatrix)

cat("\nCorrelation matrix eigenvalues:\n")
print(eigen(CorrMatrix, symmetric = TRUE)$values)


# Check whether there are missing / infinite values
cat("\nMissing values in X:\n")
print(sum(is.na(X)))

cat("\nInfinite values in X:\n")
print(sum(!is.finite(X)))

# 6. Remove SNPs with missing / infinite Z scores
keep <- apply(
  X,
  1,
  function(x) all(is.finite(x))
)

cat("\nSNPs before removing non-finite Z:", nrow(X), "\n")

X <- X[keep, , drop = FALSE]

Sumstat <- Sumstat[keep, ]

cat("SNPs after removing non-finite Z:", nrow(X), "\n")


# 7. Check correspondence between Sumstat and X
stopifnot(
  nrow(Sumstat) == nrow(X)
)

stopifnot(
  identical(
    colnames(X),
    c("Z1", "Z2")
  )
)

stopifnot(
  ncol(X) == length(Samplesize)
)


# 8. Check correlation matrix
stopifnot(
  nrow(CorrMatrix) == 2,
  ncol(CorrMatrix) == 2
)

# Make sure the order is Z1, Z2, Z3
CorrMatrix <- CorrMatrix[
  c("Z1", "Z2"),
  c("Z1", "Z2")
]

print(CorrMatrix)

# 9. Estimate Gamma distribution parameters
#    for S_Het
set.seed(12345)

cat("\nEstimating Gamma parameters...\n")

para <- EstimateGamma(
  N = 1e6,
  SampleSize = Samplesize,
  CorrMatrix = CorrMatrix
)

cat("\nGamma parameters:\n")
print(para)

names(para) <- c(
  "shape",
  "scale",
  "shift"
)

cat("\nNamed Gamma parameters:\n")
print(para)

# 10. Calculate S_Hom
cat("\nCalculating S_Hom...\n")

S_Hom <- SHom(
  X = X,
  SampleSize = Samplesize,
  CorrMatrix = CorrMatrix
)

cat("S_Hom calculation completed.\n")


# 11. Calculate P value for S_Hom
P_Hom <- pchisq(
  S_Hom,
  df = 1,
  lower.tail = FALSE
)


# 12. Calculate S_Het
#
#     Because there may be millions of SNPs,
#     calculate in blocks to reduce memory burden.

cat("\nCalculating S_Het...\n")

n <- nrow(X)

# 把全部SNP分成多少份计算
n_blocks <- 20

index_list <- split(
  seq_len(n),
  cut(
    seq_len(n),
    breaks = n_blocks,
    labels = FALSE
  )
)

S_Het <- numeric(n)

for (i in seq_along(index_list)) {
  
  cat(
    "Processing block",
    i,
    "/",
    length(index_list),
    "\n"
  )
  
  idx <- index_list[[i]]
  
  S_Het[idx] <- SHet(
    X = X[idx, , drop = FALSE],
    SampleSize = Samplesize,
    CorrMatrix = CorrMatrix
  )
}

cat("S_Het calculation completed.\n")



# 13. Calculate P value for S_Het
############################################################

P_Het <- pgamma(
  q = S_Het - para["shift"],
  shape = para["shape"],
  scale = para["scale"],
  lower.tail = FALSE
)


# 14. Check results

cat("\nSummary of S_Hom:\n")
print(summary(S_Hom))

cat("\nSummary of P_Hom:\n")
print(summary(P_Hom))

cat("\nSummary of S_Het:\n")
print(summary(S_Het))

cat("\nSummary of P_Het:\n")
print(summary(P_Het))

# 15. Construct final CPASSOC result

CPASSOC_result <- data.frame(
  
  variant_id = Sumstat$variant_id,
  
  Z1 = X[, "Z1"],
  Z2 = X[, "Z2"],
  #Z3 = X[, "Z3"],
  
  S_Hom = S_Hom,
  P_Hom = P_Hom,
  
  S_Het = S_Het,
  P_Het = P_Het
  
)

# 16. Genome-wide significant SNPs
############################################################

GWAS_sig <- CPASSOC_result[
  CPASSOC_result$P_Het < 5e-8,#不同phenotype之间遗传效应较一致用P_Hom，否则P_Het，犹豫不决用P_Het
]

cat(
  "\nNumber of SNPs with P_Het < 5e-8:",
  nrow(GWAS_sig),
  "\n"
)


# 17. Save complete CPASSOC results
CPASSOC_result[[1]] <- sub("^[^:]+:([^:]+):.*$", "\\1", CPASSOC_result[[1]])

write.csv(
  CPASSOC_result,
  "/path/to/project/CPASSOC/CPASSOC/results/CPASSOC_result.csv",
  quote = FALSE,
  row.names = FALSE
)

# 18. Save genome-wide significant results
write.table(
  GWAS_sig,
  "/path/to/project/CPASSOC/CPASSOC/results/CPASSOC_genomewide_significant.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 19. Save Gamma parameters
write.table(
  data.frame(
    parameter = names(para),
    value = as.numeric(para)
  ),
  "/path/to/project/CPASSOC/CPASSOC/results/CPASSOC_Gamma_parameters.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

