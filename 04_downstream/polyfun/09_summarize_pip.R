# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# 第 9 步：汇总 SNP 级 PIP 与可信集
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
library(dplyr)

input_dir <- "GWAS/PolyFun/finemap"
output_dir <- "GWAS/PolyFun/finemap/summary"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(input_dir, pattern = "^double_right\\.chr[0-9]+\\.[0-9]+_[0-9]+\\.gz$", full.names = TRUE)

cat("Result files:", length(files), "\n")

if (length(files) != 2748) {
  stop("Expected 2748 result files, found ", length(files))
}

# 1. Read all region results
all_data <- vector("list", length(files))

for (i in seq_along(files)) {
  cat("\rReading", i, "/", length(files), basename(files[i]))
  x <- read.delim(
    gzfile(files[i]),
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Remove any accidental malformed/header rows
  x <- x[is.numeric(x$PIP) |
           !is.na(suppressWarnings(as.numeric(x$PIP))), ]
  
  x$PIP <- as.numeric(x$PIP)
  x$P <- as.numeric(x$P)
  
  # Region name
  x$REGION <- basename(files[i])
  
  all_data[[i]] <- x
}

cat("\n")

fm <- bind_rows(all_data)

cat("Total rows:", nrow(fm), "\n")
cat("Unique SNPs:", n_distinct(fm$SNP), "\n")

# 2. SNP-level summary
snp_summary <- fm %>%
  group_by(CHR, SNP, BP, A1, A2) %>%
  summarise(
    SNPVAR = first(SNPVAR),
    N = max(N, na.rm = TRUE),
    Z = first(Z),
    P = first(P),
    PIP_MAX = max(PIP, na.rm = TRUE),
    PIP_MEAN = mean(PIP, na.rm = TRUE),
    PIP_MIN = min(PIP, na.rm = TRUE),
    BETA_MEAN = mean(BETA_MEAN, na.rm = TRUE),
    BETA_SD_MEAN = mean(BETA_SD, na.rm = TRUE),
    N_REGIONS = n(),
    N_CREDIBLE_REGIONS = sum(CREDIBLE_SET > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(PIP_MAX))

# 3. Save complete SNP-level table
write.table(snp_summary, file = gzfile(file.path(output_dir, "double_right_finemap_snp_summary.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 4. Save high-PIP SNPs
pip05 <- snp_summary %>%
  filter(PIP_MAX >= 0.5)

pip01 <- snp_summary %>%
  filter(PIP_MAX >= 0.1)

write.table(pip05, file = gzfile(file.path(output_dir, "double_right_PIP05.tsv.gz")), 
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(pip01, file = gzfile(file.path(output_dir, "double_right_PIP01.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 5. Credible-set records
credible <- fm %>%
  filter(!is.na(CREDIBLE_SET), CREDIBLE_SET > 0) %>%
  select(REGION, CHR, SNP, BP, A1, A2, P, PIP, BETA_MEAN, BETA_SD, CREDIBLE_SET) %>%
  arrange(CHR, BP, REGION, CREDIBLE_SET, desc(PIP))

write.table(credible, file = gzfile(file.path(output_dir, "double_right_credible_sets.tsv.gz")), 
            sep = "\t", quote = FALSE, row.names = F)

# 6. Summary
cat("\n========================================\n")
cat("Fine-mapping summary\n")
cat("========================================\n")

cat("Regions:             ", length(files), "\n")
cat("SNP-region rows:     ", nrow(fm), "\n")
cat("Unique SNPs:         ", nrow(snp_summary), "\n")
cat("PIP >= 0.01:         ", sum(snp_summary$PIP_MAX >= 0.01), "\n")
cat("PIP >= 0.05:         ", sum(snp_summary$PIP_MAX >= 0.05), "\n")
cat("PIP >= 0.10:         ", sum(snp_summary$PIP_MAX >= 0.10), "\n")
cat("PIP >= 0.20:         ", sum(snp_summary$PIP_MAX >= 0.20), "\n")
cat("PIP >= 0.50:         ", sum(snp_summary$PIP_MAX >= 0.50), "\n")
cat("PIP >= 0.80:         ", sum(snp_summary$PIP_MAX >= 0.80), "\n")
cat("PIP >= 0.90:         ", sum(snp_summary$PIP_MAX >= 0.90), "\n")
cat("PIP >= 0.95:         ", sum(snp_summary$PIP_MAX >= 0.95), "\n")
cat("Credible-set rows:   ", nrow(credible), "\n")

cat("\nTop 20 SNPs:\n")
print(
  snp_summary %>%
    select( CHR, SNP, BP, P,  PIP_MAX, PIP_MEAN, PIP_MIN, N_REGIONS, N_CREDIBLE_REGIONS) %>%
    head(20)
)
