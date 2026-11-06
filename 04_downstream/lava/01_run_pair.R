# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# 单个性状对的 LAVA 局部遗传相关分析（模板，改参数即可复用）
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
library(LAVA)
library(dplyr)
library(ggplot2)
library(foreach)
library(doParallel)
library(data.table)



# LAVA Hip vs Knee − Hip
input <- process.input(
  input.info.file = "/path/to/project/GWAS/LAVA/input.info.hip_double.txt",
  sample.overlap.file = "/path/to/project/GWAS/LAVA/sample.overlap.hip_double.txt",
  ref.prefix = "/path/to/ref/g1000_eur/g1000_eur",
  phenos = c("single_right_hip", "double_right")
)


loci = read.loci("/path/to/project/GWAS/LAVA/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")

locus <- process.locus(loci[1, ], input, drop.failed = FALSE)


registerDoParallel(cores = 20)
result = foreach(i = 1:nrow(loci), .packages = 'foreach') %dopar% {
  locus = process.locus(loci[i, ], input, drop.failed = F)
  
  result = list()
  if (!is.null(locus)) {
    res_univ = cbind(loci = i, run.univ(locus))
    result[[1]] = res_univ
    
    if (sum(locus$failed) == 0) {
      res_bivar = cbind(loci = i, run.bivar(locus))
      result[[2]] = res_bivar
    }
  }
  
  result
}
save(result, file = '/path/to/project/GWAS/LAVA/lava.hip.double.Rdata')