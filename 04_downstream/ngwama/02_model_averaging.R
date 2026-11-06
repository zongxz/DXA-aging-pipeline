# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# Model Averaging GWAMA（双性状）
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
# MA-GWAMA: Knee + Hip
# 2-trait Model Averaging GWAMA
library(data.table)
library(metafor)
library(AICcmodavg)
library(parallel)

# 0. File paths
knee_file <- "/path/to/project/GWAS/gwas_result/single_right_knee.txt.gz"
hip_file <- "/path/to/project/GWAS/gwas_result/single_right_hip.txt.gz"

output_dir <- "/path/to/project/GWAS/NGWAMA/result"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Read GWAS summary statistics
cat("Reading Knee GWAS...\n")
knee <- fread(knee_file)

cat("Reading Hip GWAS...\n")
hip <- fread(hip_file)

cat("\n========================================\n")
cat("Original GWAS dimensions\n")
cat("========================================\n")

cat("Knee:", nrow(knee), "rows\n")
cat("Hip :", nrow(hip), "rows\n")


# 2. Check required columns

required_cols <- c("CHR", "SNP", "POS", "A1", "A2", "N", "AF1", "BETA", "SE", "P")

missing_knee <- setdiff(required_cols, names(knee))
missing_hip  <- setdiff(required_cols, names(hip))

if (length(missing_knee) > 0) {
  stop("Missing columns in Knee GWAS: ",
       paste(missing_knee, collapse = ", "))
}

if (length(missing_hip) > 0) {
  stop("Missing columns in Hip GWAS: ",
       paste(missing_hip, collapse = ", "))
}

# 3. Remove duplicated SNPs
#
# 原始 GWAS 中存在重复 SNP ID。
#
# MA-GWAMA 后续是按照 SNP 进行 trait matching，
# 因此这里必须保证每个 GWAS 中每个 rs 只有一行。
#

cat("\n========================================\n")
cat("Removing duplicated SNPs\n")
cat("========================================\n")

cat("Knee duplicated SNP rows:", sum(duplicated(knee$SNP)), "\n")

cat("Hip duplicated SNP rows :", sum(duplicated(hip$SNP)), "\n")

knee <- knee[!duplicated(knee$SNP), ]

hip <- hip[!duplicated(hip$SNP), ]

cat("Knee after dedup:", nrow(knee), "\n")
cat("Hip  after dedup:", nrow(hip), "\n")


# 4. Standardize Knee GWAS
T1 <- data.frame(
  cptid = paste(knee$CHR, knee$POS, sep = ":"),
  rs    = knee$SNP,
  chr   = knee$CHR,
  bp    = knee$POS,
  a1    = knee$A1,
  a2    = knee$A2,
  eaf   = knee$AF1,
  n     = knee$N,
  z     = knee$BETA / knee$SE,
  p     = knee$P,
  beta  = knee$BETA,
  se    = knee$SE,
  stringsAsFactors = FALSE
)


# 5. Standardize Hip GWAS
T2 <- data.frame(
  cptid = paste(hip$CHR, hip$POS, sep = ":"),
  rs    = hip$SNP,
  chr   = hip$CHR,
  bp    = hip$POS,
  a1    = hip$A1,
  a2    = hip$A2,
  eaf   = hip$AF1,
  n     = hip$N,
  z     = hip$BETA / hip$SE,
  p     = hip$P,
  beta  = hip$BETA,
  se    = hip$SE,
  stringsAsFactors = FALSE
)

# 6. Remove non-finite values
T1 <- T1[is.finite(T1$beta) &
           is.finite(T1$se) &
           T1$se > 0 &
           is.finite(T1$eaf) &
           T1$eaf > 0 &
           T1$eaf < 1, ]

T2 <- T2[is.finite(T2$beta) &
           is.finite(T2$se) &
           T2$se > 0 &
           is.finite(T2$eaf) &
           T2$eaf > 0 &
           T2$eaf < 1, ]


# 7. Match SNPs between Knee and Hip

cat("\n========================================\n")
cat("Matching SNPs\n")
cat("========================================\n")

M1 <- merge(T1[, c("rs", "beta", "se")], T2[, c("rs", "beta", "se")], by = "rs", sort = FALSE)

# 8. Rename columns
colnames(M1) <- c("rs", "beta_knee", "se_knee", "beta_hip", "se_hip")

# 9. Remove missing values
M1 <- M1[complete.cases(M1), ]

# 10. Final duplicate check
M1 <- M1[!duplicated(M1$rs), ]

# 11. Keep SNP order consistent
M1 <- M1[match(T1$rs, M1$rs), ]
M1 <- M1[!is.na(M1$rs), ]

# 12. Number of SNPs

cat("SNPs used in MA-GWAMA:", nrow(M1), "\n")

# 13. Check allele alignment
cat("\n========================================\n")
cat("Checking allele alignment\n")
cat("========================================\n")

allele_check <- merge(T1[, c("rs", "a1", "a2")], T2[, c("rs", "a1", "a2")], by = "rs", sort = FALSE)

same_allele <- (allele_check$a1.x == allele_check$a1.y &
                  allele_check$a2.x == allele_check$a2.y)

reverse_allele <- (allele_check$a1.x == allele_check$a2.y &
                     allele_check$a2.x == allele_check$a1.y)

cat("Common SNPs:", nrow(allele_check), "\n")

cat("Same allele orientation:", sum(same_allele), "\n")

cat("Reverse allele orientation:", sum(reverse_allele), "\n")

cat("Other:", sum(!same_allele & !reverse_allele), "\n")

if (any(!same_allele)) {
  stop("Allele mismatch detected between Knee and Hip GWAS.")
  
}


# 14. Construct beta and SE matrices
#
# Rows    = SNPs
# Columns = traits
#
# Column 1 = Knee
# Column 2 = Hip
#
B <- as.matrix(M1[, c("beta_knee", "beta_hip")])
SE <- as.matrix(M1[, c("se_knee", "se_hip")])

# 15. Basic data checks
cat("\n========================================\n")
cat("Data checks\n")
cat("========================================\n")
cat("Number of SNPs:", nrow(B), "\n")
cat("Number of traits:", ncol(B), "\n")
cat("Non-finite beta:", sum(!is.finite(B)), "\n")
cat("Non-finite SE:", sum(!is.finite(SE)), "\n")
cat("SE <= 0:", sum(SE <= 0), "\n")


# 16. CTI matrix
#
# LDSC results:
#
# Knee intercept = 1.0005
# Hip  intercept = 1.0077
# Cross-trait intercept = 0.2286
#
# Trait order:
#
# 1 = Knee
# 2 = Hip
#
CTI <- matrix(c(1.0005, 0.2286, 0.2286, 1.0077),
              nrow = 2,
              byrow = TRUE)

rownames(CTI) <- c("Knee", "Hip")

colnames(CTI) <- c("Knee", "Hip")

cat("\n========================================\n")
cat("CTI matrix\n")
cat("========================================\n")

print(CTI)

# 17. SNP heritability
# LDSC:
#
# Knee h2 = 0.1017
# Hip  h2 = 0.1247
#
# MA-GWAMA official implementation uses sqrt(h2)
#

h2 <- c(sqrt(0.1017), sqrt(0.1247))

names(h2) <- c("Knee", "Hip")

cat("\n========================================\n")
cat("h2 vector used in MA-GWAMA\n")
cat("========================================\n")

print(h2)


# 18. Trait-specific model grid
# expand.grid gives:
#
#       trait1 trait2
# 1        0      0
# 2        1      0
# 3        0      1
# 4        1      1
#
# After transpose:
#
# grid[,1] = c(0,0)
# grid[,2] = c(1,0)
# grid[,3] = c(0,1)
# grid[,4] = c(1,1)
#

grid <- expand.grid(c(0, 1), c(0, 1))
grid <- t(grid)

cat("\n========================================\n")
cat("Model grid\n")
cat("========================================\n")

print(grid)


# 19. Output matrix

output_SNP <- matrix(NA_real_, nrow = nrow(B), ncol = 4)

colnames(output_SNP) <- c("Beta_Knee", "SE_Knee", "Beta_Hip", "SE_Hip")


# 20. MA-GWAMA

cat("\n========================================\n")
cat("Starting MA-GWAMA\n")
cat("========================================\n")

# Parallel MA-GWAMA

library(parallel)

# 1. Number of CPU cores

total_cores <- detectCores()

cat("Total CPU cores detected:", total_cores, "\n")

# 建议不要占满所有 CPU
# 例如服务器有 64 核，就使用 60 核
n_cores <- 32

cat("CPU cores used:", n_cores, "\n")


# 2. Define SNP-level function

MA_GWAMA_one_SNP <- function(i, B, SE, CTI, h2, grid) {
  # Sampling variance-covariance matrix
  
  V <- diag(SE[i, ]) %*%
    CTI %*%
    diag(SE[i, ])
  
  # Effect sizes
  yi <- as.numeric(B[i, ])
  
  # Model 1
  #
  # Common/shared effect
  #
  # beta_Knee = shared
  # beta_Hip  = shared
  m21 <- rma.mv(yi ~ 0 + h2, V = V, method = "ML")
  
  # Model 2
  #
  # Knee-specific effect
  #
  # beta_Knee = shared + Knee-specific
  # beta_Hip  = shared
  m22 <- rma.mv(yi ~ 0 + h2 + I(grid[, 2] * h2), V = V, method = "ML")
  
  # Candidate models
  Mods <- list(m21, m22)
  K <- c(1, 2)
  
  # Log likelihood
  LL <- lapply(Mods, logLik)
  LL <- unlist(LL)
  
  # Predictions
  XX <- sapply(Mods, predict)
  
  # Estimates
  est2 <- matrix(unlist(XX[1, 1]), 2, 1)
  est3 <- matrix(unlist(XX[1, 2]), 2, 1)
  est <- cbind(est2, est3)
  
  # SE
  se2 <- matrix(unlist(XX[2, 1]), 2, 1)
  se3 <- matrix(unlist(XX[2, 2]), 2, 1)
  se <- cbind(se2, se3)
  
  # Model averaging: Knee
  y1 <- modavgCustom(LL, K, nobs = 2, estimate = est[1, ], se = se[1, ], second.ord = TRUE)
  
  # Model averaging: Hip
  y2 <- modavgCustom(LL, K, nobs = 2, estimate = est[2, ], se = se[2, ], second.ord = TRUE)
  
  # Return
  c(y1$Mod.avg.est, y1$Uncond.SE, y2$Mod.avg.est, y2$Uncond.SE)
}

# 3. Create cluster
cat("\nCreating", n_cores, "workers...\n")
cl <- makeCluster(n_cores, type = "PSOCK")

# 4. Load required packages on workers
clusterEvalQ(cl, {
  library(metafor)
  library(AICcmodavg)
})

# 5. Export required objects
clusterExport(
  cl,
  varlist = c("B", "SE", "CTI", "h2", "grid", "MA_GWAMA_one_SNP"),
  envir = .GlobalEnv
)

# 6. Split SNPs into blocks
n_snp <- nrow(B)

# 每个 block 1000 SNP
block_size <- 1000

blocks <- split(seq_len(n_snp), ceiling(seq_len(n_snp) / block_size))
cat("Total SNPs:", n_snp, "\n")
cat("Block size:", block_size, "\n")
cat("Number of blocks:", length(blocks), "\n")

# 7. Parallel function for one block
process_block <- function(idx) {
  res <- matrix(NA_real_, nrow = length(idx), ncol = 4)
  for (j in seq_along(idx)) {
    res[j, ] <- MA_GWAMA_one_SNP(idx[j], B, SE, CTI, h2, grid)
  }
  res
}

# 8. Start parallel computation
cat("\n========================================\n")
cat("Starting parallel MA-GWAMA\n")
cat("========================================\n")

start_time <- Sys.time()
result_list <- parLapplyLB(cl, blocks, process_block)
end_time <- Sys.time()

# 9. Stop cluster
stopCluster(cl)

# 10. Combine results
output_SNP <- do.call(rbind, result_list)
colnames(output_SNP) <- c("Beta_Knee", "SE_Knee", "Beta_Hip", "SE_Hip")

# 11. Check order
if (nrow(output_SNP) != nrow(B)) {
  stop("Number of output SNPs does not match B.")
}

# 12. Runtime
cat("\n========================================\n")
cat("Parallel MA-GWAMA finished\n")
cat("========================================\n")

cat("Start time:\n")
print(start_time)

cat("End time:\n")
print(end_time)

cat("Total running time:\n")
print(end_time - start_time)

cat("Output SNPs:", nrow(output_SNP), "\n")
# 23. Save raw output_SNP
raw_output_file <- file.path(output_dir, "hip_knee_output_SNP.RData")
save(output_SNP, file = raw_output_file)

cat("Raw output saved to:\n", raw_output_file, "\n")

# 24. Match SNP annotation

A1 <- T1[match(M1$rs, T1$rs), c("rs", "cptid", "chr", "bp", "a1", "a2")]

# 25. Check annotation matching
if (nrow(A1) != nrow(output_SNP)) {
  stop("Number of SNP annotations does not match output_SNP.")
}

if (!all(A1$rs == M1$rs)) {
  stop("SNP order mismatch between annotation and MA-GWAMA results.")
}

# 26. Calculate P values
Z_Knee <- output_SNP[, "Beta_Knee"] / output_SNP[, "SE_Knee"]
Z_Hip <- output_SNP[, "Beta_Hip"] / output_SNP[, "SE_Hip"]

P_Knee <- pchisq(Z_Knee^2, df = 1, lower.tail = FALSE)
P_Hip <- pchisq(Z_Hip^2, df = 1, lower.tail = FALSE)

# 27. Construct final output
Output <- cbind(A1, output_SNP, P_Knee = P_Knee, P_Hip = P_Hip)

# 28. Rename columns
names(Output) <- c(
  "RS", "cptid", "CHR", "BP", "A1", "A2", "Beta_Knee", 
  "SE_Knee", "Beta_Hip", "SE_Hip", "P_Knee", "P_Hip"
)

# 29. Save final MA-GWAMA results
final_output_file <- file.path(output_dir, "knee_hip_MA_GWAMA.txt.gz")
fwrite(Output, final_output_file, sep = "\t")

# 30. Summary
cat("\n========================================\n")
cat("Finished\n")
cat("========================================\n")

cat("Output SNPs:", nrow(Output), "\n")

cat(
  "Genome-wide significant Knee SNPs (P < 5e-8):",
  sum(Output$P_Knee < 5e-8, na.rm = TRUE),
  "\n"
)

cat(
  "Genome-wide significant Hip SNPs (P < 5e-8):",
  sum(Output$P_Hip < 5e-8, na.rm = TRUE),
  "\n"
)

cat("Final output:\n", final_output_file, "\n")

cat("\nFirst 10 results:\n")

print(head(Output, 10))