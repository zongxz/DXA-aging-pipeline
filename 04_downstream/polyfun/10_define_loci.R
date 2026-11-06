# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# 第 10 步：按 250 kb 断点划分 locus 并汇总
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
library(dplyr)

input <- "GWAS/PolyFun/finemap/summary/double_right_finemap_snp_summary.tsv.gz"
outdir <- "GWAS/PolyFun/finemap/summary"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

fm <- read.delim(
  gzfile(input),
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 1. Keep PIP >= 0.05
x <- fm %>%
  filter(PIP_MAX >= 0.05) %>%
  arrange(CHR, BP)

cat("PIP >= 0.05 SNPs:", nrow(x), "\n")

# 2. Define loci
#
# New locus if:
#   chromosome changes
#   OR distance from previous SNP > 250 kb

x <- x %>%
  mutate(
    prev_chr = lag(CHR),
    prev_bp = lag(BP),
    
    new_locus = ifelse(is.na(prev_chr) |
                         CHR != prev_chr |
                         (BP - prev_bp) > 250000, 1, 0),
    
    LOCUS = cumsum(new_locus)
  )

# 3. Locus summary
locus_summary <- x %>%
  group_by(LOCUS, CHR) %>%
  summarise(
    START = min(BP),
    END = max(BP),
    
    N_SNP = n(),
    
    LEAD_SNP = SNP[which.max(PIP_MAX)],
    LEAD_BP = BP[which.max(PIP_MAX)],
    
    PIP_MAX = max(PIP_MAX),
    
    LEAD_P = P[which.max(PIP_MAX)],
    
    N_CREDIBLE_SNP = sum(N_CREDIBLE_REGIONS > 0, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(CHR, START)

# 4. Save
write.table(
  locus_summary,
  file = gzfile(file.path(
    outdir, "double_right_finemap_loci.tsv.gz"
  )),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  x,
  file = gzfile(
    file.path(outdir, "double_right_PIP05_with_locus.tsv.gz")
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 5. Print summary
cat("\n========================================\n")
cat("Fine-mapped loci\n")
cat("========================================\n")

cat("PIP >= 0.05 SNPs:", nrow(x), "\n")
cat("Loci:", nrow(locus_summary), "\n")

cat("\nLocus summary:\n")

print(locus_summary %>%
        arrange(desc(PIP_MAX)))
